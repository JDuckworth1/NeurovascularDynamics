clc;clear;close all;
DataManager = FileManager;
dataset_name = 'DKLab';
stack_list = {'20240129R152'}; %'20231219R155Ec'
num_stack = numel(stack_list);
dev_mode = true;
%%    
if ~dev_mode
    pool_obj = gcp('nocreate');
    delete(pool_obj);
    num_processes = 8;
    parpool(num_processes);
    num_thread_per_prop = 48 / num_processes;
end
%%
for stack_idx = 1 : num_stack
    %%
    stack = stack_list{stack_idx};
    
    info_fp = fullfile(DataManager.fp_raw_data_folder(dataset_name, stack), ...
        sprintf('%s_%s_raw_data_info.mat', dataset_name, stack));
    
    data_info = DataManager.load_data(info_fp);       
    vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Method');
    %% Parameters
    para = struct;
    para.linescan_rate_Hz = data_info.parameters.linescan_rate_Hz;
    para.pixel_size_um = data_info.parameters.pixel_size_um;
    % reshape the array - make frame rate consistent
    para.target_frame_rate_Hz = 50;
    para.num_lines_per_frame = round(para.linescan_rate_Hz / para.target_frame_rate_Hz);
    para.frame_rate_Hz = para.linescan_rate_Hz / para.num_lines_per_frame;
    para.num_init_line_dropped = 20;
    %
    para.vis_Q = false && dev_mode;
    para.vis_bg_est_Q = true && dev_mode;
    submedfiltQ = 1; %ADDED 1/17/24 to correct for GFP
        
    % Global background estimation
    para.bg.est_bg_int_n = 0.2;
    para.bg.est_bg_snr = 3;
    
    % Vessel image segmentation
    para.vm.disk_r = 3;
    para.vm.sph_r = 3;
    para.vm.min_cc_size = 10;
    
    % Radius estimation
    para.re.min_snr = 3;
    
    % Flux estimation
    % Fraction of pixels in the vessel that are background in the RBC
    % channel
    para.fe.est_bg_frac = 0.9;
    para.fe.min_rbc_size_pxl = 10;
    %% 
    %parpool(10)
    parfor iter_file = 1 : data_info.num_files
        roi_idx = data_info.roi_list(iter_file);
        record_idx = data_info.record_idx_list(iter_file);
        if ~dev_mode
            maxNumCompThreads(num_thread_per_prop);
        end
        tic_roi = tic;
        
        result_str = struct;
        result_str.dataset_name = dataset_name;
        result_str.stack = stack;
        result_str.roi_idx = roi_idx;
        result_str.record_idx = record_idx;
        
        result_str.data_fn_base = sprintf('%s_%s_roi_%03d_%05d', dataset_name, ...
            stack, roi_idx, record_idx);
        
        result_str.radius_data_fn = sprintf('%s_vsl.mat', result_str.data_fn_base);
        result_str.rbc_data_fn = sprintf('%s_rbc.mat', result_str.data_fn_base);
        result_str.radius_data_fp = fullfile(DataManager.fp_raw_data_folder(dataset_name, stack), ...
            result_str.radius_data_fn);
        result_str.flow_data_fp = fullfile(DataManager.fp_raw_data_folder(dataset_name, stack), ...
            result_str.rbc_data_fn);
        
        %Modified 1/27/24   
        result_str.result_fp = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack), ...
            sprintf('%s_%s_roi_%03d_%05d_sf.mat', dataset_name, stack,...
            roi_idx, record_idx));
        %% Load data
        try
            vsl_data = DataManager.load_data(result_str.radius_data_fp);
            flow_data = DataManager.load_data(result_str.flow_data_fp);
            vsl_im = permute(vsl_data.im, [2,1,3]);
            ls_im = permute(flow_data.im, [2,1,3]);
        catch ME
            fprintf('Failed to load data in ROI %d recording %d. Skip...\n', ...
                roi_idx, record_idx);
            if ~dev_mode
                %continue;
            end
        end
        %%
        try
            %% Set up parameters
            para_str = para;
            % Reshape frames
            raw_im_sz = size(vsl_im);
            para_str.num_frame = floor((raw_im_sz(2) * raw_im_sz(3) -  para_str.num_init_line_dropped)...
                / para_str.num_lines_per_frame);
            para_str.num_lines = para_str.num_frame * para_str.num_lines_per_frame;
            para_str.num_pixel_per_line = raw_im_sz(1);
            
            vsl_im = reshape(vsl_im, raw_im_sz(1), []);
            vsl_im = vsl_im(:, (para_str.num_init_line_dropped + 1) : ...
                (para_str.num_init_line_dropped + para_str.num_lines));
            vsl_im = reshape(vsl_im, raw_im_sz(1), para_str.num_lines_per_frame, []);

            
            ls_im = reshape(ls_im, raw_im_sz(1), []);
            ls_im = ls_im(:, (para_str.num_init_line_dropped + 1) : ...
                (para_str.num_init_line_dropped + para_str.num_lines));
            ls_im = reshape(ls_im, raw_im_sz(1), para_str.num_lines_per_frame, []);
            
            if strcmp(stack,'JD230627M997_101123') && ismember(roi_idx,[16,17,18,19,20]) %Crop this frame b/c we picked up a dural vessel on the edge
                vsl_im = vsl_im(1:50,:,:); 
                ls_im = ls_im(1:50,:,:); 
            end
            
            %Try median-correcting as estimate of noise in lumen.
            if submedfiltQ == 1
                bkg_med = squeeze(median(ls_im,2));
                bkg_med = repmat(bkg_med,[1,1,316]);
                bkg_med = permute(bkg_med,[1,3,2]);
                ls_im = ls_im - bkg_med;
            end
            
            para_str.im_size = [para_str.num_pixel_per_line, para_str.num_lines_per_frame];            
            para_str.dt_s = 1 ./ para_str.frame_rate_Hz; 
            %% Smooth the image
            for i = 1 : para_str.num_frame
                vsl_im(:, :, i) = medfilt2(vsl_im(:, :, i)); %Default
                %vsl_im(:, :, i) = medfilt2(vsl_im(:, :, i),[5,5]);
               % ls_im(:, :, i) = medfilt2(ls_im(:, :, i)); 
               %ls_im(:,:,i) = imgaussfilt(ls_im(:,:,i));
            end
            disp('done smoothing')
            %% Global background estiamtion
            tmp_bg_tic = tic;
            est_bg_str = fun_DKLab_TBRL_background_estiamtion(vsl_im, para_str.bg.est_bg_int_n, ...
                para_str.bg.est_bg_snr, para_str.vis_bg_est_Q);
            fprintf('Finish estimating background level. Elapsed time is %.2f seconds\n', ...
                toc(tmp_bg_tic));            
            %% Compute vessel mask 
            vsl_mask = VRRBCFlux.compute_vessel_mask(vsl_im, est_bg_str.est_bg_int, ...
                para_str.vm.min_cc_size, para_str.vm.disk_r, para_str.vm.sph_r);
            %% RBC segmentation 
            rbc_seg_cell = VRRBCFlux.compute_rbc_cc_ReaChR(ls_im, vsl_mask, ...
                para_str.fe.est_bg_frac, para_str.fe.min_rbc_size_pxl);            
            rbc_prop = VRRBCFlux.analyze_single_frame_rbc_cc(rbc_seg_cell, ...
                ls_im);
            %%  Radius estimation
            r_str = VRRBCFlux.estimate_radius_for_each_frame(vsl_im, vsl_mask, ... 
                para_str.re.min_snr, para_str.pixel_size_um);            
            [line_r_trace, line_r_str] = VRRBCFlux.estimate_radius_for_each_line(vsl_im, vsl_mask, ... 
                para_str.re.min_snr, para_str.pixel_size_um);
            %% Save results
            result_str.parameters = para_str;
            result_str.global_background = est_bg_str;
            result_str.radius = r_str;
            result_str.radius_l_um = line_r_trace;
%             result_str.radius_l_str = line_r_str; % This file is too
%             large...
            result_str.flux = rbc_prop;
            
            DataManager.write_data(result_str.result_fp, result_str);
        catch ME
            fprintf('Failed to process ROI %d record %d. Error message: %s\n', ...
                roi_idx, record_idx, ME.message);
            if ~dev_mode
                continue
            end
        end
        fprintf('Finish processing ROI %d record %d. Elapsed time is %.2f seconds\n', ...
            roi_idx, record_idx, toc(tic_roi));
    end
end
%%
% med_rbc_trace = reshape(rbc_prop.rbc_trace, para_str.num_lines_per_frame, []);
% raw_rbc_trace = reshape(rbc_prop_raw.rbc_trace, para_str.num_lines_per_frame, []);
% med_rbc_trace_sm = movmean(sum(med_rbc_trace, 1), 50);
% raw_rbc_trace_sm = movmean(sum(raw_rbc_trace, 1), 50);
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% plot(ax_hdl, raw_rbc_trace_sm);
% hold(ax_hdl, 'on');
% plot(ax_hdl, med_rbc_trace_sm);
%%
rbc_mask = cellfun(@(x) x.rbc_mask, rbc_seg_cell, 'UniformOutput', false);
rbc_mask = cat(3, rbc_mask{:});
implay(rbc_mask);
%%
vis_idx = 100;
fig_hdl = figure;
a1 = nexttile();
imagesc(a1, ls_im(:, :, vis_idx));
a1.DataAspectRatio = [1,1,1];
a2 = nexttile();
imagesc(a2, rbc_mask(:, :, vis_idx));
a2.DataAspectRatio = [1,1,1];
a3 = nexttile();
imshowpair(ls_im(:, :, vis_idx), rbc_mask(:, :, vis_idx));
%% 
figure;
start = round(1200);
endpt = start+100;
for i = start:endpt
   subplot(2,1,1)
   imagesc(ls_im(:,:,i))
   caxis([0 2.5*10^4]);
   subplot(2,1,2)
   imagesc(rbc_seg_cell{i,1}.rbc_mask)
%     imagesc(vsl_im(:,:,i))
   title(sprintf('%.0f',i));
   pause()
end

%%


