clc;clear;close all;
DataManager = FileManager;
dataset_name = 'DKLab';
%stack_list = {'TBRLBB20211129', 'TBRLBB20211119', 'TBRLBB20211118',...
    %'TBRLBB20211103', 'TBRLBB20211108'}
%stack_list = {'JD230627M002_92723','JD230627M002_100423','JD230627M997_100123','JD230627M997_100723','JD230627M997_101123','JD230627M945_101123'};
stack_list = {'20231219R155Ec','JD230523R1495trials_011524','JD230523R149_011524','20240124R152','20240129R155','20240129R152'};

% int_stack_name = 'TBRLBB';
int_stack_name = 'JDRLXJ';
int_fig_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    int_stack_name), 'Integrated_sl');
num_stack = numel(stack_list);
dev_mode = false;
vis_Q = false;
% Processing parameters
ppp_str = struct;
ppp_str.num_phase_bin = 16;
ppp_str.num_flux_bin = 16;
ppp_str.num_r_bin = 16;
ppp_str.num_speed_bin = 16;
ppp_str.phase_bin_edge = linspace(-pi, pi, ppp_str.num_phase_bin + 1);
ppp_str.vasomotion_min_fraction = 0.25;
pps_str.frac_psd_f_range_Hz = [5, 12];
%% Load all the data 
stack_data = cell(num_stack, 1);
for stack_idx = 1 : num_stack
    %%
    stack = stack_list{stack_idx};
    data_info = DataManager.load_data(fullfile(DataManager.fp_raw_data_folder(dataset_name, ...
        stack), sprintf('%s_%s_raw_data_info.mat', dataset_name, stack)));
    switch stack %See Lab book p64
        case '20231219R155Ec'
            roi_idx_pia = [4, 5, 7, 8];
            roi_idx_ptr = [];
            roi_idx_piaV = [];
            vis_data_idx = [95, 95, 95, 95];
        case 'JD230523R1495trials_011524'
            roi_idx_pia = [10, 14, 15];
            roi_idx_ptr = [];
            roi_idx_piaV = [];
            vis_data_idx = [97, 97, 97];            
        case 'JD230523R149_011524' 
            roi_idx_pia = [8, 9, 13];
            roi_idx_ptr = [];
            roi_idx_piaV = [];
            vis_data_idx = [95, 97, 97];            
        case '20240124R152' 
            roi_idx_pia = [3, 4, 5, 6, 7, 8, 13];
            roi_idx_ptr = [9, 10, 11];
            roi_idx_piaV = [];
            vis_data_idx = [98, 98, 98, 98, 98, 98, 98];            
        case '20240129R155'
            roi_idx_pia = [11, 14, 16];
            roi_idx_ptr = [];
            roi_idx_piaV = [];
            vis_data_idx = [98, 98, 98];
        case '20240129R152'
            roi_idx_pia = [18, 20];
            roi_idx_ptr = [];
            roi_idx_piaV = [];
            vis_data_idx = [98, 99];
    end
    [roi_list_idx, roi_idx] = fun_bin_data_to_idx_list(data_info.roi_list);
    num_roi = numel(roi_idx);
    tmp_roi_data = cell(num_roi, 1);
    exp_data = cell(length(vis_data_idx), 1);

    % Load data
%     for iter_file = 1 : data_info.num_files
    for iter_file = 1:length(vis_data_idx)
        tic_vis = tic;
        roi_idx = roi_idx_pia(iter_file);
        %roi_idx = data_info.roi_list(iter_file);
        %record_idx = data_info.record_idx_list(iter_file);
        record_idx = 1; %Use first trial unless otherwise noted!
        vis_num = vis_data_idx(iter_file);
        
        fprintf('Processing ROI %d record %d\n', roi_idx, record_idx);
        result_fp = fullfile(DataManager.fp_analysis_visdata_folder(dataset_name, stack, vis_num), ...
            'vis_data', sprintf('%s_%s_ROI_%03d_%05d_vis_data.mat',...
            dataset_name, stack, roi_idx, record_idx));
        try
            vis_data = DataManager.load_data(result_fp);
        catch ME
            warning('Failed to load data for ROI %d', roi_idx);
            vis_data = [];
            continue;
        end
        exp_data{iter_file} = vis_data;
    end

    [~, roi_idx] = fun_bin_data_to_idx_list(data_info.roi_list);
    num_roi = numel(roi_idx);
    tmp_roi_data = cell(num_roi, 1);
    roi_list_idx = 1:length(vis_data_idx);
        
      %% Plot R2 vs slope histograms
    for iter_roi = 1:length(vis_data_idx)
        tmp_roi_idx = roi_idx_pia(iter_roi);
        %tmp_roi_idx = roi_idx(iter_roi); %ROI number
        if any(roi_idx_pia == tmp_roi_idx)
            tmp_vs_type = 'Pia';
        elseif any(roi_idx_ptr == tmp_roi_idx)
            tmp_vs_type = 'Penetrating';
        elseif any(roi_idx_piaV == tmp_roi_idx)
            tmp_vs_type = 'PiaV';
        else %If this trial is excluding
            tmp_vs_type = [];
        end            
        
        tmp_list_idx = roi_list_idx(iter_roi);
        tmp_data = cat(1, exp_data{tmp_list_idx});
        tmp_data_cell = exp_data(tmp_list_idx(1):tmp_list_idx(end));
      
        tmp_hdl = VRRFAI(tmp_data, tmp_vs_type);        
        tmp_hdl.parse_parameters(ppp_str);
        tmp_hdl.merge_radius_data();
        tmp_hdl.merge_flux_and_speed_data();
        
        %% Bin data according to hphase 
        r_name = 'r_or_mm';
        flux_name = 'per_frame_mm';
        
        if ~contains(cell2mat(stack_list(stack_idx)),'5trials')
            r_trace = tmp_hdl.radius.(r_name).data;
            r_stat = tmp_hdl.radius.(r_name).stat;
            r_trace_n = tmp_hdl.radius.(r_name).data_n;
            r_edge = linspace(r_stat.prctile_val(2), r_stat.prctile_val(end-1),...
                tmp_hdl.parameter.num_r_edge);
            
            flux_trace = tmp_hdl.flux.(flux_name).data;
            flux_stat = tmp_hdl.flux.(flux_name).stat;
            flux_trace_n = tmp_hdl.flux.(flux_name).data_n;
            flux_stat_n = tmp_hdl.flux.(flux_name).stat_n;
            flux_edge = linspace(flux_stat.prctile_val(2), flux_stat.prctile_val(end-1),...
                tmp_hdl.parameter.num_flux_edge);
            flux_n_edge = linspace(flux_stat_n.prctile_val(2), flux_stat_n.prctile_val(end-1),...
                tmp_hdl.parameter.num_flux_edge);
        else %For 5 trials, just use first trial for this analysis.
            tmp_data_cell = exp_data(tmp_list_idx(1):tmp_list_idx(end));
            
            %Get individual radius/flux traces
            if strcmp(tmp_vs_type,'Pia') 
            tmp_data_ind = cellfun(@(x) VRRFAI(x,'Pia'), tmp_data_cell, 'UniformOutput',false);
            else
            tmp_data_ind = cellfun(@(x) VRRFAI(x,'PiaV'), tmp_data_cell, 'UniformOutput',false);   
            end
            tmp_hdl = tmp_data_ind{1};
            tmp_hdl.parse_parameters(ppp_str);
            tmp_hdl.merge_radius_data();
            tmp_hdl.merge_flux_and_speed_data();
%             r_trace_ind = cell(1,length(tmp_list_idx));
%             flux_trace_ind = cell(1,length(tmp_list_idx));
%             tmp_data_ind{1}.parse_parameters(ppp_str);
%             tmp_data_ind{1}.merge_radius_data();
%             tmp_data_ind{1}.merge_flux_and_speed_data();
            r_trace = tmp_hdl.radius.('r_or_mm').data; %Just take first trial from these
            flux_trace = tmp_hdl.flux.('per_frame_mm').data; %Just take first trial
            %tmp_hdl = tmp_data_ind{1};
            
            
            r_trace = tmp_hdl.radius.(r_name).data;
            r_stat = tmp_hdl.radius.(r_name).stat;
            r_trace_n = tmp_hdl.radius.(r_name).data_n;
            r_edge = linspace(r_stat.prctile_val(2), r_stat.prctile_val(end-1),...
                tmp_hdl.parameter.num_r_edge);
            
            flux_trace = tmp_hdl.flux.(flux_name).data;
            flux_stat = tmp_hdl.flux.(flux_name).stat;
            flux_trace_n = tmp_hdl.flux.(flux_name).data_n;
            flux_stat_n = tmp_hdl.flux.(flux_name).stat_n;
            flux_edge = linspace(flux_stat.prctile_val(2), flux_stat.prctile_val(end-1),...
                tmp_hdl.parameter.num_flux_edge);
            flux_n_edge = linspace(flux_stat_n.prctile_val(2), flux_stat_n.prctile_val(end-1),...
                tmp_hdl.parameter.num_flux_edge);
            
        end
        
%         phase_trace = tmp_hdl.radius.r_or_mm.phase_l;
%         phase_edge = tmp_hdl.parameter.phase_bin_edge;
        
        speed_trace = tmp_hdl.speed.(flux_name).data;
        speed_stat = tmp_hdl.speed.(flux_name).stat;
        speed_trace_n = tmp_hdl.speed.(flux_name).data_n;
        speed_edge = linspace(speed_stat.min, speed_stat.max, ...
            tmp_hdl.parameter.num_flux_edge);
        speed_n_stat = tmp_hdl.speed.(flux_name).stat_n;
        speed_n_edge = linspace(speed_n_stat.prctile_val(2), speed_n_stat.prctile_val(end-1), ...
            tmp_hdl.parameter.num_flux_edge);
        
%         vflux_trace = tmp_hdl.speed.(flux_name).est_flux;
%         vflux_trace_n = tmp_hdl.speed.(flux_name).est_flux_n;
%         vflux_stat = fun_analysis_get_basic_statistics(tmp_hdl.speed.(flux_name).est_flux);
%         vflux_edge = linspace(vflux_stat.prctile_val(2), vflux_stat.prctile_val(end-1), ...
%             tmp_hdl.parameter.num_flux_edge);
%         vflux_n_stat = fun_analysis_get_basic_statistics(vflux_trace_n);
%         vflux_n_edge = linspace(vflux_n_stat.prctile_val(2), vflux_n_stat.prctile_val(end-1), ...
%             tmp_hdl.parameter.num_flux_edge);  
        
        %Calculate restricted range here for traces!! then use them below.
        lowp = designfilt('lowpassiir','PassbandFrequency',0.5,'StopbandFrequency',0.75,'SampleRate',50,'PassbandRipple',4,'StopbandAttenuation',10);
        r_tracefilt = filtfilt(lowp,double(r_trace));
%         figure; subplot(2,1,1); plot(r_trace); hold on; plot(r_tracefilt);
%         subplot(2,1,2); plot(diff(diff(r_tracefilt)));
        d2r = diff(diff(r_tracefilt));
        [d2max,d2maxloc] = max(d2r(1:20*50)); %Look between 0 and 20 seconds
        r_trace_infl{1} = r_trace(1:(d2maxloc+1)); %Use these like r_trace is used.
        r_trace_infl{2} = r_trace((d2maxloc+2):end);
        r_trace_infl_n{1} = (r_trace_infl{1}-mean(r_trace_infl{1}))./mean(r_trace_infl{1});
        r_trace_infl_n{2} = (r_trace_infl{2}-mean(r_trace_infl{2}))./mean(r_trace_infl{2});
        flux_trace_infl{1} = flux_trace(1:(d2maxloc+1));
        flux_trace_infl{2} = flux_trace((d2maxloc+2):end);
        flux_trace_infl_n{1} = (flux_trace_infl{1}-mean(flux_trace_infl{1}))./mean(flux_trace_infl{1});
        flux_trace_infl_n{2} = (flux_trace_infl{2}-mean(flux_trace_infl{2}))./mean(flux_trace_infl{2});        
        %Use constant time as well (shouldn't change results)
        tlim = 15;
        r_trace_t{1} = r_trace(1:(50*tlim)); %Use these like r_trace is used.
        r_trace_t{2} = r_trace((50*tlim+1):end);
        r_trace_t_n{1} = (r_trace_t{1}-mean(r_trace_t{1}))./mean(r_trace_t{1});
        r_trace_t_n{2} = (r_trace_t{2}-mean(r_trace_t{2}))./mean(r_trace_t{2});        
        flux_trace_t{1} = flux_trace(1:(50*tlim));
        flux_trace_t{2} = flux_trace((50*tlim+1):end);         
        flux_trace_t_n{1} = (flux_trace_t{1}-mean(flux_trace_t{1}))./mean(flux_trace_t{1});
        flux_trace_t_n{2} = (flux_trace_t{2}-mean(flux_trace_t{2}))./mean(flux_trace_t{2});  
        %Use 20pct decrease as separation point.
        pctlim = 0.20;
        r_20pct = r_trace(1) - pctlim*r_trace(1);
        r_diff = r_trace - r_20pct;
        rmin_vals = find(r_diff < 0); %Use first passing
        if ~isempty(rmin_vals)
            r20pt = rmin_vals(1);
            r_trace_20pt{1} = r_trace(1:r20pt); %Use these like r_trace is used.
            r_trace_20pt{2} = r_trace((r20pt+1):end);
            r_trace_20pt_n{1} = (r_trace_20pt{1}-mean(r_trace_20pt{1}))./mean(r_trace_20pt{1});
            r_trace_20pt_n{2} = (r_trace_20pt{2}-mean(r_trace_20pt{2}))./mean(r_trace_20pt{2});
            flux_trace_20pt{1} = flux_trace(1:r20pt);
            flux_trace_20pt{2} = flux_trace((r20pt+1):end);
            flux_trace_20pt_n{1} = (flux_trace_20pt{1}-mean(flux_trace_20pt{1}))./mean(flux_trace_20pt{1});
            flux_trace_20pt_n{2} = (flux_trace_20pt{2}-mean(flux_trace_20pt{2}))./mean(flux_trace_20pt{2});
            %TraceCheck_JD230523R149_011524_1_24_24
%             pathstr = ''
%             cd(pathstr)
%             figure;
%             t = (0:length(r_trace)-1)/50+1/50/2;
%             plot(t,r_trace); ylabel('Radius (um)','Interpreter','latex');
%             yline(r_20pct)
%             hold on; plot(t(1:r20pt),r_trace(1:r20pt),'r');
%             yyaxis right
%             plot(t,flux_trace); ylabel('Flux (/s)','Interpreter','latex');
%             title(tmp_hdl.title_prefix,'Interpreter','latex')
%             xlabel('time (s)','Interpreter','latex');
%             savefig([tmp_hdl.title_prefix,'_20pctLine'])
%             saveas(gcf,[tmp_hdl.title_prefix,'_20pctLine.png']);
            
          tmp_hdl.radius.('r_or_mm').stat.mean = mean(r_trace_20pt{1});
        end
        
        %% Radius vs flux - original 
        tmp_hdl.radius_vs_flux = struct;
        tmp_hdl.radius_vs_flux.raw = tmp_hdl.analyze_x_vs_y(r_trace,...
            flux_trace);
        if vis_Q
            fig_hdl = tmp_hdl.vis_radius_vs_flux(r_trace, flux_trace, ...
                r_edge, flux_edge, tmp_hdl.radius_vs_flux.raw, 'Radius (\mum)', ...
                'Flux (RBC/s)', tmp_hdl.title_prefix);
            tmp_hdl.save_figure(fig_hdl, {}, 'radius_vs_flux.png', ~dev_mode);
        end
        %% Radius vs flux - normalized
        tmp_hdl.radius_vs_flux.normalized = tmp_hdl.analyze_normalized_changes(...
            r_trace_n, flux_trace_n);   
        tmp_hdl.radius_vs_flux.normalized_intercept = tmp_hdl.analyze_normalized_changes_intercept(...
            r_trace_n, flux_trace_n);         
        if vis_Q
            fig_hdl = tmp_hdl.vis_radius_vs_flux(r_trace_n, flux_trace_n, ...
                tmp_hdl.radius_vs_flux.normalized.x_bin_edge,...
                flux_n_edge, tmp_hdl.radius_vs_flux.normalized, '\Deltar/\langler\rangle', ...
                '\Deltaq/\langleq\rangle', tmp_hdl.title_prefix);
            tmp_hdl.save_figure(fig_hdl, {}, 'dr_vs_dq_n.png', ~dev_mode);
        end
        %% Radius vs flux - normalized - smaller range
        %inflection point
        tmp_hdl.radius_vs_flux.norm_infl1 = tmp_hdl.analyze_normalized_changes(...
            r_trace_infl_n{1},flux_trace_infl_n{1});
        tmp_hdl.radius_vs_flux.norm_infl2 = tmp_hdl.analyze_normalized_changes(...
            r_trace_infl_n{2},flux_trace_infl_n{2});
        %constant time
        tmp_hdl.radius_vs_flux.norm_t1 = tmp_hdl.analyze_normalized_changes(...
            r_trace_t_n{1},flux_trace_t_n{1});
        tmp_hdl.radius_vs_flux.norm_t2 = tmp_hdl.analyze_normalized_changes(...
            r_trace_t_n{2},flux_trace_t_n{2});   
        %20pct drop
        if ~isempty(rmin_vals)
            %Save norm traces
            tmp_hdl.norm_r = r_trace_20pt_n{1};
            tmp_hdl.norm_flux = flux_trace_20pt_n{1};
            tmp_hdl.radius_vs_flux.norm_20pt1 = tmp_hdl.analyze_normalized_changes(...
                r_trace_20pt_n{1},flux_trace_20pt_n{1}); %Initial segment
            tmp_hdl.radius_vs_flux.norm_20pt2 = tmp_hdl.analyze_normalized_changes(...
                r_trace_20pt_n{2},flux_trace_20pt_n{2}); %Final segment
            tmp_hdl.radius_vs_flux.norm_20pt1_intercept = tmp_hdl.analyze_normalized_changes_intercept(...
                r_trace_20pt_n{1},flux_trace_20pt_n{1});   %Initial segment             
        end
        if vis_Q
            fig_hdl = tmp_hdl.vis_radius_vs_flux(r_trace_20pt_n{1}, flux_trace_20pt_n{1}, ...
                tmp_hdl.radius_vs_flux.norm_20pt1.x_bin_edge,...
                flux_n_edge, tmp_hdl.radius_vs_flux.norm_20pt1, '\Deltar/\langler\rangle', ...
                '\Deltaq/\langleq\rangle', tmp_hdl.title_prefix);
            tmp_hdl.save_figure(fig_hdl, {}, 'dr_vs_dq_n_infl.png', ~dev_mode);
        end
        if vis_Q
            fig_hdl = tmp_hdl.vis_radius_vs_flux(r_trace_t_n{1}, flux_trace_t_n{1}, ...
                tmp_hdl.radius_vs_flux.norm_t1.x_bin_edge,...
                flux_n_edge, tmp_hdl.radius_vs_flux.norm_t1, '\Deltar/\langler\rangle', ...
                '\Deltaq/\langleq\rangle', tmp_hdl.title_prefix);
            tmp_hdl.save_figure(fig_hdl, {}, 'dr_vs_dq_n_infl.png', ~dev_mode);
        end        
        
        
        
        %% Radius vs speed - raw
        tmp_hdl.radius_vs_speed = struct;
        tmp_hdl.radius_vs_speed.raw = tmp_hdl.analyze_x_vs_y(r_trace,...
            speed_trace);
        if vis_Q
            fig_hdl = tmp_hdl.vis_radius_vs_flux(r_trace, speed_trace, ...
                r_edge, speed_edge, tmp_hdl.radius_vs_speed.raw, 'Radius (\mum)', ...
                'Speed (mm/s)', tmp_hdl.title_prefix);
            tmp_hdl.save_figure(fig_hdl, {}, 'radius_vs_speed.png', ~dev_mode);
        end
        %% Radius vs speed - normalized
        tmp_hdl.radius_vs_speed.normalized = tmp_hdl.analyze_normalized_changes(r_trace_n,...
            speed_trace_n);
        if vis_Q
            fig_hdl = tmp_hdl.vis_radius_vs_flux(r_trace_n, speed_trace_n, ...
                tmp_hdl.radius_vs_speed.normalized.x_bin_edge,...
                speed_n_edge, tmp_hdl.radius_vs_speed.normalized, '\Deltar/\langler\rangle', ...
                '\Deltav/\langlev\rangle', tmp_hdl.title_prefix);
            tmp_hdl.save_figure(fig_hdl, {}, 'radius_vs_speed_n.png', ~dev_mode);
        end
        %% Flux vs volumetric flow - normalized
%         tmp_hdl.flux_vs_flow = struct;
%         tmp_hdl.flux_vs_flow.normalized = tmp_hdl.analyze_normalized_changes(vflux_trace_n,...
%             flux_trace_n);
        
        %% Flux vs volumetric flow - original
%         tmp_hdl.flux_vs_flow.raw = tmp_hdl.analyze_x_vs_y(vflux_trace, flux_trace);        

%         %% Analyze vasodilation intervals - flux vs radius
%         % Outlier rejections
%         vm_int_stat_flux = tmp_hdl.select_vasomition_interval_stat(...
%             tmp_hdl.flux.(flux_name).vm_int_stat, ...
%             tmp_hdl.parameter.vasomotion_min_fraction);
% %         if vis_Q
% %             fig_hdl = tmp_hdl.vis_vasomotion_interval_stat_v2(vm_int_stat_flux,...
% %                 sprintf('%s %s', tmp_hdl.title_prefix, 'interval flux vs radius'));
% %             tmp_hdl.save_figure(fig_hdl, {}, 'vasomotion_r2q_fit_stat.png', ~dev_mode);
% %         end
%         %% Analyze vasodilation intervals - speed vs radius
%         % Outlier rejections
%         vm_int_stat_speed = tmp_hdl.select_vasomition_interval_stat(...
%             tmp_hdl.speed.(flux_name).vm_int_stat,...
%             tmp_hdl.parameter.vasomotion_min_fraction);
% %         if vis_Q
% %             fig_hdl = tmp_hdl.vis_vasomotion_interval_stat_v2(vm_int_stat_speed,...
% %                 sprintf('%s %s', tmp_hdl.title_prefix, 'interval speed vs radius'));
% %             tmp_hdl.save_figure(fig_hdl, {}, 'vasomotion_r2v_fit_stat.png', ~dev_mode);
% %         end
        %% Frequency
        % Heart beat frequency @ 50 Hz resolution 
%         tmp_hdl.frequency = struct;
%         tmp_hdl.frequency.fpsd_gt_5Hz = tmp_hdl.compute_fractional_power_in_frequency_range(pps_str.frac_psd_f_range_Hz);        
        % High pass filter and then 
        drop_time_points = 100;        
        tmp_hdl.compute_radius_modulation_in_frequency_range(pps_str.frac_psd_f_range_Hz, 100);
        %%

        tmp_roi_data{iter_roi} = tmp_hdl;
        if isempty(rmin_vals)
            tmp_roi_data{iter_roi} = [];
        end
    end

    %%
    stack_data{stack_idx} = tmp_roi_data;
end
%% Integration
stack_data = cat(1, stack_data{:});
stack_data = cat(1, stack_data{:});

vsl_type_list = arrayfun(@(x) x.vessel_type, stack_data, 'UniformOutput', false);
is_pia_vessel_Q = strcmpi(vsl_type_list, 'Pia');
is_piaV_vessel_Q = strcmpi(vsl_type_list, 'PiaV');
is_PA_vessel_Q = ~is_pia_vessel_Q & ~is_piaV_vessel_Q;
vessel_avg_r_list = arrayfun(@(x) x.radius.(r_name).stat.mean, stack_data);

all_norm_r = arrayfun(@(x) x.norm_r, stack_data, 'UniformOutput', false);
all_norm_r = cell2mat(all_norm_r);
all_norm_flux = arrayfun(@(x) x.norm_flux, stack_data, 'UniformOutput', false);
all_norm_flux = cell2mat(all_norm_flux);
%% 
% is_pia_vessel_Q([19,25]) = [];
% is_PA_vessel_Q([19,25]) = [];
% #########################################################################
%% David's figures
fit_fieldname = {'fit_raw', 'fit_med'};
fit_title_name = {'Data', 'Binned histogram median'};
num_fn = numel(fit_fieldname);
    %% Radius vs flux 
    for i_p = 1 : num_fn
        tmp_fieldname = fit_fieldname{i_p};
        tmp_title = fit_title_name{i_p};
        q2r_k_list = arrayfun(@(x) x.radius_vs_flux.normalized.(tmp_fieldname).k, stack_data);
        q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.normalized.(tmp_fieldname).R2, stack_data);

        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        sc_hdl = scatter(ax_hdl, q2r_k_list(is_pia_vessel_Q), q2r_fit_r2_list(is_pia_vessel_Q), ...
            150, 2*vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');
%         hold(ax_hdl, 'on');
%         sc_hdl_1 = scatter(ax_hdl, q2r_k_list(is_PA_vessel_Q), q2r_fit_r2_list(is_PA_vessel_Q), ...
%             150, 2*vessel_avg_r_list(is_PA_vessel_Q), '*');
%         sc_hdl_2 = scatter(ax_hdl, q2r_k_list(is_piaV_vessel_Q), q2r_fit_r2_list(is_piaV_vessel_Q), ...
%             150, 2*vessel_avg_r_list(is_piaV_vessel_Q), 'square');   %ADDED FOR VEINS, ReaChR EXPERIMENTS   

        stat_str = struct;
        %Pia stats
        stat_str.pia_a2c_a = mean(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_std = std(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wa = sum(q2r_k_list(is_pia_vessel_Q) .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c2_wa = sum(q2r_k_list(is_pia_vessel_Q).^2 .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wstd = sqrt(stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2);
        df_fac = sum(q2r_fit_r2_list(is_pia_vessel_Q))^2/((sum(q2r_fit_r2_list(is_pia_vessel_Q))^2) ...
            - sum(q2r_fit_r2_list(is_pia_vessel_Q).^2));
        stat_str.pia_a2c_wstd_corrected = sqrt((stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2)*df_fac);
        %PA stats
%         stat_str.pa_a2c_a = mean(q2r_k_list(is_PA_vessel_Q));
%         stat_str.pa_a2c_std = std(q2r_k_list(is_PA_vessel_Q));    
%         stat_str.pa_a2c_wa = sum(q2r_k_list(is_PA_vessel_Q) .* q2r_fit_r2_list(is_PA_vessel_Q)) / ...
%             sum(q2r_fit_r2_list(is_PA_vessel_Q));
%         stat_str.pa_a2c2_wa = sum(q2r_k_list(is_PA_vessel_Q).^2 .* q2r_fit_r2_list(is_PA_vessel_Q)) / ...
%             sum(q2r_fit_r2_list(is_PA_vessel_Q));
%         stat_str.pa_a2c_wstd = sqrt(stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2); 
%         df_fac = sum(q2r_fit_r2_list(is_PA_vessel_Q))^2/((sum(q2r_fit_r2_list(is_PA_vessel_Q))^2) ...
%             - sum(q2r_fit_r2_list(is_PA_vessel_Q).^2));
%         stat_str.pa_a2c_wstd_corrected = sqrt((stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2)*df_fac);
%         %Vein stats
%         stat_str.piaV_a2c_a = mean(q2r_k_list(is_piaV_vessel_Q));
%         stat_str.piaV_a2c_std = std(q2r_k_list(is_piaV_vessel_Q));    
%         stat_str.piaV_a2c_wa = sum(q2r_k_list(is_piaV_vessel_Q) .* q2r_fit_r2_list(is_piaV_vessel_Q)) / ...
%             sum(q2r_fit_r2_list(is_piaV_vessel_Q));
%         stat_str.piaV_a2c2_wa = sum(q2r_k_list(is_piaV_vessel_Q).^2 .* q2r_fit_r2_list(is_piaV_vessel_Q)) / ...
%             sum(q2r_fit_r2_list(is_piaV_vessel_Q));
%         stat_str.piaV_a2c_wstd = sqrt(stat_str.piaV_a2c2_wa - stat_str.piaV_a2c_wa^2); 
%         df_fac = sum(q2r_fit_r2_list(is_piaV_vessel_Q))^2/((sum(q2r_fit_r2_list(is_piaV_vessel_Q))^2) ...
%             - sum(q2r_fit_r2_list(is_piaV_vessel_Q).^2));
%         stat_str.piaV_a2c_wstd_corrected = sqrt((stat_str.piaV_a2c2_wa - stat_str.piaV_a2c_wa^2)*df_fac);

%         legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
%             stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd_corrected), ...
%             sprintf('PA\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
%             stat_str.pa_a2c_a, stat_str.pa_a2c_std, stat_str.pa_a2c_wa, stat_str.pa_a2c_wstd_corrected), ...
%             sprintf('PiaV\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
%             stat_str.piaV_a2c_a, stat_str.piaV_a2c_std, stat_str.piaV_a2c_wa, stat_str.piaV_a2c_wstd_corrected)};
        legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd_corrected)};
        legend(ax_hdl, [sc_hdl], legend_str, 'Location', 'best');

        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Normalized flux to radius slope';
        ax_hdl.YLabel.String = 'Fitting R^2';
        ax_hdl.YLim = [0,1];
        sc_hdl.MarkerFaceAlpha = 0.6;
        %sc_hdl_1.MarkerFaceAlpha = 0.6;
        c_bar = colorbar(ax_hdl);
        cm_hdl = colormap(ax_hdl, 'jet');
        %c_bar.Label.String = 'Radius (\mum)';
        c_bar.Label.String = 'Diameter (\mum)';
        ax_hdl.FontSize = 12;
        ax_hdl.Title.String = tmp_title;

        fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_flux_vs_radius_fit_vsD%s_stat_ReaChR_2_4_24.png', ...
            dataset_name, int_stack_name, tmp_fieldname));
%         fun_print_image_in_several_formats(fig_hdl, fig_fp);

        data_table = struct;
        data_table.q2r_k = q2r_k_list;
        data_table.q2r_fit_R2 = q2r_fit_r2_list;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
%         writetable(data_table, table_fp);
    end
%     tmpmat = [q2r_k_list,q2r_fit_r2_list];
%     for i = 1:size(tmpmat,1)
%        tmpmat(i,3) = str2num(cell2mat(extractBetween(exp_data{i}.vis_fp_template,'ROI_','_00001')));
%     end
        %% Radius vs flux, initial decay
   
        for i_p = 2:2
            tmp_fieldname = fit_fieldname{i_p};
            tmp_title = fit_title_name{i_p};
            q2r_k_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1.(tmp_fieldname).k, stack_data);
            q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1.(tmp_fieldname).R2, stack_data);
            
            
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            sc_hdl = scatter(ax_hdl, q2r_k_list(is_pia_vessel_Q), q2r_fit_r2_list(is_pia_vessel_Q), ...
                150, 2*vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');
            %         hold(ax_hdl, 'on');
            %         sc_hdl_1 = scatter(ax_hdl, q2r_k_list(is_PA_vessel_Q), q2r_fit_r2_list(is_PA_vessel_Q), ...
            %             150, 2*vessel_avg_r_list(is_PA_vessel_Q), '*');
            %         sc_hdl_2 = scatter(ax_hdl, q2r_k_list(is_piaV_vessel_Q), q2r_fit_r2_list(is_piaV_vessel_Q), ...
            %             150, 2*vessel_avg_r_list(is_piaV_vessel_Q), 'square');
            %             %ADDED FOR VEINS, ReaChR EXPERIMENTS, don't plot these for
            %             now.
            
            stat_str = struct;
            %Pia stats
            stat_str.pia_a2c_a = mean(q2r_k_list(is_pia_vessel_Q));
            stat_str.pia_a2c_std = std(q2r_k_list(is_pia_vessel_Q));
            stat_str.pia_a2c_wa = sum(q2r_k_list(is_pia_vessel_Q) .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
                sum(q2r_fit_r2_list(is_pia_vessel_Q));
            stat_str.pia_a2c2_wa = sum(q2r_k_list(is_pia_vessel_Q).^2 .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
                sum(q2r_fit_r2_list(is_pia_vessel_Q));
            stat_str.pia_a2c_wstd = sqrt(stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2);
            df_fac = sum(q2r_fit_r2_list(is_pia_vessel_Q))^2/((sum(q2r_fit_r2_list(is_pia_vessel_Q))^2) ...
                - sum(q2r_fit_r2_list(is_pia_vessel_Q).^2));
            stat_str.pia_a2c_wstd_corrected = sqrt((stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2)*df_fac);
            %PA stats
            %         stat_str.pa_a2c_a = mean(q2r_k_list(is_PA_vessel_Q));
            %         stat_str.pa_a2c_std = std(q2r_k_list(is_PA_vessel_Q));
            %         stat_str.pa_a2c_wa = sum(q2r_k_list(is_PA_vessel_Q) .* q2r_fit_r2_list(is_PA_vessel_Q)) / ...
            %             sum(q2r_fit_r2_list(is_PA_vessel_Q));
            %         stat_str.pa_a2c2_wa = sum(q2r_k_list(is_PA_vessel_Q).^2 .* q2r_fit_r2_list(is_PA_vessel_Q)) / ...
            %             sum(q2r_fit_r2_list(is_PA_vessel_Q));
            %         stat_str.pa_a2c_wstd = sqrt(stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2);
            %         df_fac = sum(q2r_fit_r2_list(is_PA_vessel_Q))^2/((sum(q2r_fit_r2_list(is_PA_vessel_Q))^2) ...
            %             - sum(q2r_fit_r2_list(is_PA_vessel_Q).^2));
            %         stat_str.pa_a2c_wstd_corrected = sqrt((stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2)*df_fac);
            %         %Vein stats
            %         stat_str.piaV_a2c_a = mean(q2r_k_list(is_piaV_vessel_Q));
            %         stat_str.piaV_a2c_std = std(q2r_k_list(is_piaV_vessel_Q));
            %         stat_str.piaV_a2c_wa = sum(q2r_k_list(is_piaV_vessel_Q) .* q2r_fit_r2_list(is_piaV_vessel_Q)) / ...
            %             sum(q2r_fit_r2_list(is_piaV_vessel_Q));
            %         stat_str.piaV_a2c2_wa = sum(q2r_k_list(is_piaV_vessel_Q).^2 .* q2r_fit_r2_list(is_piaV_vessel_Q)) / ...
            %             sum(q2r_fit_r2_list(is_piaV_vessel_Q));
            %         stat_str.piaV_a2c_wstd = sqrt(stat_str.piaV_a2c2_wa - stat_str.piaV_a2c_wa^2);
            %         df_fac = sum(q2r_fit_r2_list(is_piaV_vessel_Q))^2/((sum(q2r_fit_r2_list(is_piaV_vessel_Q))^2) ...
            %             - sum(q2r_fit_r2_list(is_piaV_vessel_Q).^2));
            %         stat_str.piaV_a2c_wstd_corrected = sqrt((stat_str.piaV_a2c2_wa - stat_str.piaV_a2c_wa^2)*df_fac);
            
            %         legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            %             stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd_corrected), ...
            %             sprintf('PA\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            %             stat_str.pa_a2c_a, stat_str.pa_a2c_std, stat_str.pa_a2c_wa, stat_str.pa_a2c_wstd_corrected), ...
            %             sprintf('PiaV\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            %             stat_str.piaV_a2c_a, stat_str.piaV_a2c_std, stat_str.piaV_a2c_wa, stat_str.piaV_a2c_wstd_corrected)};
            legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
                stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd_corrected)};
            legend(ax_hdl, [sc_hdl], legend_str, 'Location', 'best');
            
            grid(ax_hdl, 'on');
            ax_hdl.XLabel.String = 'Normalized flux to radius slope';
            ax_hdl.YLabel.String = 'Fitting R^2';
            ax_hdl.YLim = [0,1];
            sc_hdl.MarkerFaceAlpha = 0.6;
            %         sc_hdl_1.MarkerFaceAlpha = 0.6;
            c_bar = colorbar(ax_hdl);
            cm_hdl = colormap(ax_hdl, 'jet');
            %c_bar.Label.String = 'Radius (\mum)';        fun_print_image_in_several_formats(fig_hdl, fig_fp);
            
            c_bar.Label.String = 'Diameter (\mum)';
            ax_hdl.FontSize = 12;
            ax_hdl.Title.String = tmp_title;
            
            xlim([0 4]);
            pbaspect([1.42,1,1]);
            xline(stat_str.pia_a2c_wa)
            
            fig_fp = fullfile(int_fig_folder, [sprintf('%s_%s_flux_vs_radius_fit_vsD%s_stat_ReaChR_InitialSeg_20pt_2_4_24_DaspectUpdate_Xline.png', ...
                dataset_name, int_stack_name, tmp_fieldname)]);
            %fig_fp = fullfile(int_fig_folder, ['20240124R152',sprintf('%s_%s_flux_vs_radius_fit_vsD%s_stat_ReaChR_InitialSeg_20pt_1_25_24.png', ...
            %dataset_name, int_stack_name, tmp_fieldname)]);
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            
            data_table = struct;
            data_table.q2r_k = q2r_k_list;
            data_table.q2r_fit_R2 = q2r_fit_r2_list;
            data_table.avg_radius_um = vessel_avg_r_list;
            data_table.is_pia = is_pia_vessel_Q;
            data_table = struct2table(data_table);
            table_fp = strrep(fig_fp, 'png', 'csv');tmp_fieldname = fit_fieldname{i_p};
            tmp_title = fit_title_name{i_p};
            q2r_k_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1.(tmp_fieldname).k, stack_data);
            q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1.(tmp_fieldname).R2, stack_data);
            writetable(data_table, table_fp);
        end

%% Radius vs flux, after initial drop.
   
    for i_p = 1 : num_fn
        tmp_fieldname = fit_fieldname{i_p};
        tmp_title = fit_title_name{i_p};
        q2r_k_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt2.(tmp_fieldname).k, stack_data);
        q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt2.(tmp_fieldname).R2, stack_data);

        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        sc_hdl = scatter(ax_hdl, q2r_k_list(is_pia_vessel_Q), q2r_fit_r2_list(is_pia_vessel_Q), ...
            150, 2*vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');
        hold(ax_hdl, 'on');
        sc_hdl_1 = scatter(ax_hdl, q2r_k_list(is_PA_vessel_Q), q2r_fit_r2_list(is_PA_vessel_Q), ...
            150, 2*vessel_avg_r_list(is_PA_vessel_Q), '*');
%         sc_hdl_2 = scatter(ax_hdl, q2r_k_list(is_piaV_vessel_Q), q2r_fit_r2_list(is_piaV_vessel_Q), ...
%             150, 2*vessel_avg_r_list(is_piaV_vessel_Q), 'square');

        stat_str = struct;
        %Pia stats
        stat_str.pia_a2c_a = mean(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_std = std(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wa = sum(q2r_k_list(is_pia_vessel_Q) .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c2_wa = sum(q2r_k_list(is_pia_vessel_Q).^2 .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wstd = sqrt(stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2);
        df_fac = sum(q2r_fit_r2_list(is_pia_vessel_Q))^2/((sum(q2r_fit_r2_list(is_pia_vessel_Q))^2) ...
            - sum(q2r_fit_r2_list(is_pia_vessel_Q).^2));
        stat_str.pia_a2c_wstd_corrected = sqrt((stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2)*df_fac);
        %PA stats
        stat_str.pa_a2c_a = mean(q2r_k_list(is_PA_vessel_Q));
        stat_str.pa_a2c_std = std(q2r_k_list(is_PA_vessel_Q));    
        stat_str.pa_a2c_wa = sum(q2r_k_list(is_PA_vessel_Q) .* q2r_fit_r2_list(is_PA_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_PA_vessel_Q));
        stat_str.pa_a2c2_wa = sum(q2r_k_list(is_PA_vessel_Q).^2 .* q2r_fit_r2_list(is_PA_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_PA_vessel_Q));
        stat_str.pa_a2c_wstd = sqrt(stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2); 
        df_fac = sum(q2r_fit_r2_list(is_PA_vessel_Q))^2/((sum(q2r_fit_r2_list(is_PA_vessel_Q))^2) ...
            - sum(q2r_fit_r2_list(is_PA_vessel_Q).^2));
        stat_str.pa_a2c_wstd_corrected = sqrt((stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2)*df_fac);
        %Vein stats
        stat_str.piaV_a2c_a = mean(q2r_k_list(is_piaV_vessel_Q));
        stat_str.piaV_a2c_std = std(q2r_k_list(is_piaV_vessel_Q));    
        stat_str.piaV_a2c_wa = sum(q2r_k_list(is_piaV_vessel_Q) .* q2r_fit_r2_list(is_piaV_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_piaV_vessel_Q));
        stat_str.piaV_a2c2_wa = sum(q2r_k_list(is_piaV_vessel_Q).^2 .* q2r_fit_r2_list(is_piaV_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_piaV_vessel_Q));
        stat_str.piaV_a2c_wstd = sqrt(stat_str.piaV_a2c2_wa - stat_str.piaV_a2c_wa^2); 
        df_fac = sum(q2r_fit_r2_list(is_piaV_vessel_Q))^2/((sum(q2r_fit_r2_list(is_piaV_vessel_Q))^2) ...
            - sum(q2r_fit_r2_list(is_piaV_vessel_Q).^2));
        stat_str.piaV_a2c_wstd_corrected = sqrt((stat_str.piaV_a2c2_wa - stat_str.piaV_a2c_wa^2)*df_fac);

        legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd_corrected), ...
            sprintf('PA\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pa_a2c_a, stat_str.pa_a2c_std, stat_str.pa_a2c_wa, stat_str.pa_a2c_wstd_corrected), ...
            sprintf('PiaV\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.piaV_a2c_a, stat_str.piaV_a2c_std, stat_str.piaV_a2c_wa, stat_str.piaV_a2c_wstd_corrected)};
        legend(ax_hdl, [sc_hdl,sc_hdl_1], legend_str, 'Location', 'best');

        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Normalized flux to radius slope, final segment';
        ax_hdl.YLabel.String = 'Fitting R^2';
        ax_hdl.YLim = [0,1];
        sc_hdl.MarkerFaceAlpha = 0.6;
        sc_hdl_1.MarkerFaceAlpha = 0.6;
        c_bar = colorbar(ax_hdl);
        cm_hdl = colormap(ax_hdl, 'jet');
        %c_bar.Label.String = 'Radius (\mum)';
        c_bar.Label.String = 'Diameter (\mum)';
        ax_hdl.FontSize = 12;
        ax_hdl.Title.String = tmp_title;

        fig_fp = fullfile(int_fig_folder, ['30Arterioles',sprintf('%s_%s_flux_vs_radius_fit_vsD%s_stat_ReaChR_FinalSeg_20pt_1_25_24.png', ...
            dataset_name, int_stack_name, tmp_fieldname)]);
        fun_print_image_in_several_formats(fig_hdl, fig_fp);

        data_table = struct;
        data_table.q2r_k = q2r_k_list;
        data_table.q2r_fit_R2 = q2r_fit_r2_list;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
        writetable(data_table, table_fp);
    end        
    
%% Radius vs flux with intercept, plot slope
    for i_p = 1 : num_fn
        tmp_fieldname = fit_fieldname{i_p};
        tmp_title = fit_title_name{i_p};
        q2r_k_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).k, stack_data);
        q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).R2, stack_data);
        q2r_k_SE_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).k_SE, stack_data);
        
        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        sc_hdl = scatter(ax_hdl, q2r_k_list(is_pia_vessel_Q), q2r_fit_r2_list(is_pia_vessel_Q), ...
            150, 2*vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');

        stat_str = struct;
        stat_str.pia_a2c_a = mean(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_std = std(q2r_k_list(is_pia_vessel_Q));

        stat_str.pia_a2c_wa = sum(q2r_k_list(is_pia_vessel_Q) .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c2_wa = sum(q2r_k_list(is_pia_vessel_Q).^2 .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wstd = sqrt(stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2);

        stat_str.pa_a2c_a = mean(q2r_k_list(~is_pia_vessel_Q));
        stat_str.pa_a2c_std = std(q2r_k_list(~is_pia_vessel_Q));    
        stat_str.pa_a2c_wa = sum(q2r_k_list(~is_pia_vessel_Q) .* q2r_fit_r2_list(~is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(~is_pia_vessel_Q));
        stat_str.pa_a2c2_wa = sum(q2r_k_list(~is_pia_vessel_Q).^2 .* q2r_fit_r2_list(~is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(~is_pia_vessel_Q));
        stat_str.pa_a2c_wstd = sqrt(stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2);    

        legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd)};
        legend(ax_hdl, [sc_hdl], legend_str, 'Location', 'best');

        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Normalized flux to radius slope';
        ax_hdl.YLabel.String = 'Fitting R^2';
        ax_hdl.YLim = [0,1];
        sc_hdl.MarkerFaceAlpha = 0.6;
        sc_hdl_1.MarkerFaceAlpha = 0.6;
        c_bar = colorbar(ax_hdl);
        cm_hdl = colormap(ax_hdl, 'jet');
        %c_bar.Label.String = 'Radius (\mum)';
        c_bar.Label.String = 'Diameter (\mum)';
        ax_hdl.FontSize = 12;
        ax_hdl.Title.String = tmp_title;
        xlim([0 4]);

        int_fig_folder = '/net/birdstore/Jacob/Flux/DKLab/JDRLXJ/visualization/Q_D_Fit_Intercepts';
        fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_flux_vs_radius_fit_slope_vsD_ReaChR%s_stat.png', ...
            dataset_name, int_stack_name, tmp_fieldname));
        %fun_print_image_in_several_formats(fig_hdl, fig_fp);

        data_table = struct;
        data_table.q2r_k = q2r_k_list;
        data_table.q2r_fit_R2 = q2r_fit_r2_list;
        data_table.q2r_k_SE = q2r_k_SE_list;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
        writetable(data_table, table_fp);
    end    
%% Radius vs flux with intercept, plot intercept
    for i_p = 1 : num_fn
        tmp_fieldname = fit_fieldname{i_p};
        tmp_title = fit_title_name{i_p};
        q2r_k_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).int, stack_data);
        q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).R2, stack_data);
        q2r_b_SE_list = arrayfun(@(x) x.radius_vs_flux.normalized_intercept.(tmp_fieldname).int_SE, stack_data);
        
        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        sc_hdl = scatter(ax_hdl, q2r_k_list(is_pia_vessel_Q), q2r_fit_r2_list(is_pia_vessel_Q), ...
            150, 2*vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');
        hold(ax_hdl, 'on');

        stat_str = struct;
        stat_str.pia_a2c_a = mean(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_std = std(q2r_k_list(is_pia_vessel_Q));

        stat_str.pia_a2c_wa = sum(q2r_k_list(is_pia_vessel_Q) .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c2_wa = sum(q2r_k_list(is_pia_vessel_Q).^2 .* q2r_fit_r2_list(is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wstd = sqrt(stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2);

        stat_str.pa_a2c_a = mean(q2r_k_list(~is_pia_vessel_Q));
        stat_str.pa_a2c_std = std(q2r_k_list(~is_pia_vessel_Q));    
        stat_str.pa_a2c_wa = sum(q2r_k_list(~is_pia_vessel_Q) .* q2r_fit_r2_list(~is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(~is_pia_vessel_Q));
        stat_str.pa_a2c2_wa = sum(q2r_k_list(~is_pia_vessel_Q).^2 .* q2r_fit_r2_list(~is_pia_vessel_Q)) / ...
            sum(q2r_fit_r2_list(~is_pia_vessel_Q));
        stat_str.pa_a2c_wstd = sqrt(stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2);    

        legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd), ...
            sprintf('PA\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pa_a2c_a, stat_str.pa_a2c_std, stat_str.pa_a2c_wa, stat_str.pa_a2c_wstd)};
        legend(ax_hdl, [sc_hdl], legend_str, 'Location', 'best');

        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Normalized flux to radius slope';
        ax_hdl.YLabel.String = 'Fitting R^2';
        ax_hdl.YLim = [0,1];
        sc_hdl.MarkerFaceAlpha = 0.6;
        sc_hdl_1.MarkerFaceAlpha = 0.6;
        c_bar = colorbar(ax_hdl);
        cm_hdl = colormap(ax_hdl, 'jet');
        %c_bar.Label.String = 'Radius (\mum)';
        c_bar.Label.String = 'Diameter (\mum)';
        ax_hdl.FontSize = 12;
        ax_hdl.Title.String = tmp_title;

        int_fig_folder = '/net/birdstore/Jacob/Flux/DKLab/JDRLXJ/visualization/Q_D_Fit_Intercepts';
        fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_flux_vs_radius_fit_intercept_vsD_ReaChR%s_stat.png', ...
            dataset_name, int_stack_name, tmp_fieldname));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);

        data_table = struct;
        data_table.q2r_b = q2r_k_list;
        data_table.q2r_fit_R2 = q2r_fit_r2_list;
        data_table.q2r_b_SE = q2r_b_SE_list;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
        writetable(data_table, table_fp);
    end         
    
%% Radius vs flux with intercept, plot slope vs inverse variance
    for i_p = 1 : num_fn
        tmp_fieldname = fit_fieldname{i_p};
        tmp_title = fit_title_name{i_p};
        q2r_k_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).k, stack_data);
        q2r_fit_r2_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).R2, stack_data);
        q2r_k_SE_list = arrayfun(@(x) x.radius_vs_flux.norm_20pt1_intercept.(tmp_fieldname).k_SE, stack_data);
        q2r_invVAR_list = (1./q2r_k_SE_list).^2;
        
        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        sc_hdl = scatter(ax_hdl, q2r_k_list(is_pia_vessel_Q), q2r_invVAR_list(is_pia_vessel_Q), ...
            150, 2*vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');
        hold(ax_hdl, 'on');
        sc_hdl_1 = scatter(ax_hdl, q2r_k_list(~is_pia_vessel_Q), q2r_invVAR_list(~is_pia_vessel_Q), ...
            150, 2*vessel_avg_r_list(~is_pia_vessel_Q), '*');

        stat_str = struct;
        stat_str.pia_a2c_a = mean(q2r_k_list(is_pia_vessel_Q));
        stat_str.pia_a2c_std = std(q2r_k_list(is_pia_vessel_Q));
        npia = sum(is_pia_vessel_Q);
        npa = sum(~is_pia_vessel_Q);

        stat_str.pia_a2c_wa = sum(q2r_k_list(is_pia_vessel_Q) .* q2r_invVAR_list(is_pia_vessel_Q)) / ...
            sum(q2r_invVAR_list(is_pia_vessel_Q));
        stat_str.pia_a2c2_wa = sum(q2r_k_list(is_pia_vessel_Q).^2 .* q2r_invVAR_list(is_pia_vessel_Q)) / ...
            sum(q2r_invVAR_list(is_pia_vessel_Q));
        stat_str.pia_a2c_wstd = sqrt((stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2)*npia/(npia-1));

        stat_str.pa_a2c_a = mean(q2r_k_list(~is_pia_vessel_Q));
        stat_str.pa_a2c_std = std(q2r_k_list(~is_pia_vessel_Q));    
        stat_str.pa_a2c_wa = sum(q2r_k_list(~is_pia_vessel_Q) .* q2r_invVAR_list(~is_pia_vessel_Q)) / ...
            sum(q2r_invVAR_list(~is_pia_vessel_Q));
        stat_str.pa_a2c2_wa = sum(q2r_k_list(~is_pia_vessel_Q).^2 .* q2r_invVAR_list(~is_pia_vessel_Q)) / ...
            sum(q2r_invVAR_list(~is_pia_vessel_Q));
        stat_str.pa_a2c_wstd = sqrt((stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2)*npa/(npa-1));    

        legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd), ...
            sprintf('PA\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pa_a2c_a, stat_str.pa_a2c_std, stat_str.pa_a2c_wa, stat_str.pa_a2c_wstd)};
        legend(ax_hdl, [sc_hdl, sc_hdl_1], legend_str, 'Location', 'best');

        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Normalized flux to radius slope';
        ax_hdl.YLabel.String = 'Inverse Variance';
%         ax_hdl.YLim = [0,1];
        sc_hdl.MarkerFaceAlpha = 0.6;
        sc_hdl_1.MarkerFaceAlpha = 0.6;
        c_bar = colorbar(ax_hdl);
        cm_hdl = colormap(ax_hdl, 'jet');
        %c_bar.Label.String = 'Radius (\mum)';
        c_bar.Label.String = 'Diameter (\mum)';
        ax_hdl.FontSize = 12;
        ax_hdl.Title.String = tmp_title;
        xlim([0 4])

        int_fig_folder = '/net/birdstore/Jacob/Flux/DKLab/JDRLXJ/visualization/Q_D_Fit_Intercepts';
        fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_flux_vs_radius_fit_slope_vsD_invVAR_ReaChR%s_stat.png', ...
            dataset_name, int_stack_name, tmp_fieldname));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);

        data_table = struct;
        data_table.q2r_k = q2r_k_list;
        data_table.q2r_fit_R2 = q2r_fit_r2_list;
        data_table.q2r_k_SE = q2r_k_SE_list;
        data_table.q2r_invVAR = q2r_invVAR_list;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
        writetable(data_table, table_fp);
    end     
        
    
    
    
    %% Flux to volumetric flow
    for i_p = 1 : num_fn
        tmp_fieldname = fit_fieldname{i_p};
        tmp_title = fit_title_name{i_p};
        q2Q_k_list = arrayfun(@(x) x.flux_vs_flow.normalized.(tmp_fieldname).k, stack_data);
        q2Q_fit_r2_list = arrayfun(@(x) x.flux_vs_flow.normalized.(tmp_fieldname).R2, stack_data);

        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        sc_hdl = scatter(ax_hdl, q2Q_k_list(is_pia_vessel_Q), q2Q_fit_r2_list(is_pia_vessel_Q), ...
            150, vessel_avg_r_list(is_pia_vessel_Q), 'o', 'filled');
        hold(ax_hdl, 'on');
        sc_hdl_1 = scatter(ax_hdl, q2Q_k_list(is_PA_vessel_Q), q2Q_fit_r2_list(is_PA_vessel_Q), ...
            150, vessel_avg_r_list(is_PA_vessel_Q), '*');

        stat_str = struct;
        stat_str.pia_a2c_a = mean(a2c_list(is_pia_vessel_Q));
        stat_str.pia_a2c_std = std(a2c_list(is_pia_vessel_Q));

        stat_str.pia_a2c_wa = sum(a2c_list(is_pia_vessel_Q) .* a2c_fit_r2(is_pia_vessel_Q)) / ...
            sum(a2c_fit_r2(is_pia_vessel_Q));
        stat_str.pia_a2c2_wa = sum(a2c_list(is_pia_vessel_Q).^2 .* a2c_fit_r2(is_pia_vessel_Q)) / ...
            sum(a2c_fit_r2(is_pia_vessel_Q));
        stat_str.pia_a2c_wstd = sqrt(stat_str.pia_a2c2_wa - stat_str.pia_a2c_wa^2);

        stat_str.pa_a2c_a = mean(a2c_list(is_PA_vessel_Q));
        stat_str.pa_a2c_std = std(a2c_list(is_PA_vessel_Q));    
        stat_str.pa_a2c_wa = sum(a2c_list(is_PA_vessel_Q) .* a2c_fit_r2(is_PA_vessel_Q)) / ...
            sum(a2c_fit_r2(is_PA_vessel_Q));
        stat_str.pa_a2c2_wa = sum(a2c_list(is_PA_vessel_Q).^2 .* a2c_fit_r2(is_PA_vessel_Q)) / ...
            sum(a2c_fit_r2(is_PA_vessel_Q));
        stat_str.pa_a2c_wstd = sqrt(stat_str.pa_a2c2_wa - stat_str.pa_a2c_wa^2);    

        legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pia_a2c_a, stat_str.pia_a2c_std, stat_str.pia_a2c_wa, stat_str.pia_a2c_wstd), ...
            sprintf('PA\n\t%.3f \\pm %.3f\n\t%.3f \\pm %.3f', ...
            stat_str.pa_a2c_a, stat_str.pa_a2c_std, stat_str.pa_a2c_wa, stat_str.pa_a2c_wstd)};
        legend(ax_hdl, [sc_hdl, sc_hdl_1], legend_str, 'Location', 'best');    

        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Normalized flux to volumetric flux slope';
        ax_hdl.YLabel.String = 'Fitting R^2';
        ax_hdl.YLim = [0,1];
        sc_hdl.MarkerFaceAlpha = 0.6;
        sc_hdl_1.MarkerFaceAlpha = 0.6;
        c_bar = colorbar(ax_hdl);
        cm_hdl = colormap(ax_hdl, 'jet');
        c_bar.Label.String = 'Radius (\mum)';
        ax_hdl.FontSize = 12;
        ax_hdl.Title.String = tmp_title;

        fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_flux_vs_volumetric_flux_fit_%s_stat.png', ...
            dataset_name, int_stack_name, tmp_fieldname));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);

        data_table = struct;
        data_table.q2Q_k = q2Q_k_list;
        data_table.q2Q_fit_R2 = q2Q_fit_r2_list;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
        writetable(data_table, table_fp);
    end
%% Plot ReaChR 2D histogram across all vessels
%20pt drop
x_stat = fun_analysis_get_basic_statistics(all_norm_r);
y_stat = fun_analysis_get_basic_statistics(all_norm_flux);
x_edge = linspace(x_stat.prctile_val(2), x_stat.prctile_val(end-1),17);
y_edge = linspace(y_stat.prctile_val(2), y_stat.prctile_val(end-1),17);
dyn_in_dxn_bin_str = fun_analysis_get_y_stat_in_x_bin(all_norm_r, all_norm_flux, x_edge);
tmp_fit_y = dyn_in_dxn_bin_str.y_median; %Median of y in every x bin

tmp_fit_x = dyn_in_dxn_bin_str.x_bin_val; %X bin values
tmp_fit = fitlm(tmp_fit_x, tmp_fit_y, 'Intercept', false);
dyn_in_dxn_bin_str.fit_med.k = tmp_fit.Coefficients.Estimate(1);
dyn_in_dxn_bin_str.fit_med.R2 = tmp_fit.Rsquared.Adjusted;

tmp_fit_all = fitlm(all_norm_r, all_norm_flux, 'Intercept', false);
dyn_in_dxn_bin_str.fit_raw.k = tmp_fit_all.Coefficients.Estimate(1);
dyn_in_dxn_bin_str.fit_raw.R2 = tmp_fit_all.Rsquared.Adjusted;

%% PLOT! From VRRFAI vis_xy_joint_dist_in_hist
fig_hdl = figure;
ax_hdl = axes(fig_hdl);

hist_hdl = histogram2(ax_hdl, all_norm_r, all_norm_flux, 'XBinedges', x_edge,'YBinEdges',y_edge,'DisplayStyle', 'tile');
%hist_hdl = histogram2(ax_hdl, all_norm_r, all_norm_flux,'DisplayStyle', 'tile');
hold(ax_hdl, 'on');
binned_stat = dyn_in_dxn_bin_str;
eb_hdl = errorbar(ax_hdl, binned_stat.x_bin_val, binned_stat.y_median, ...
    binned_stat.y_median - binned_stat.y_prctile(:, 3), ...
    binned_stat.y_prctile(:, 5) - binned_stat.y_median, ...
    'LineWidth', 1.5);
ax_hdl.XLabel.String = 'dD/D';
ax_hdl.YLabel.String = 'dQ/Q';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.XLabel.Interpreter = 'latex';
xlim([-0.15 0.2])
ylim([-0.5 0.6])
grid off;
%Plot fit
xfit = x_edge(1):0.01:x_edge(end);
yfit = xfit.*dyn_in_dxn_bin_str.fit_med.k;
hold on
plot(xfit,yfit,'k','LineWidth',2)
ax_hdl.TickLabelInterpreter = 'latex';

str = sprintf('k = %.2f, R$^2$ = %.2f',dyn_in_dxn_bin_str.fit_med.k,dyn_in_dxn_bin_str.fit_med.R2);
title({'ReaChR experiments, initial 20pct drop joint distribution',str},'Interpreter','latex');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
caxis([0 200])
%%
cd('');
print(fig_hdl,'-depsc','-painters','ReaChR_JointDist_AllTrials');

    %% Heart rate frequency contribution
    plt_fn = {'r', 'flux', 'speed'};
    for i_p = 1 : numel(plt_fn)
        tmp_fn = plt_fn{i_p};
        tmp_fpsd_r = arrayfun(@(x) x.frequency.fpsd_gt_5Hz.(tmp_fn).stat.median, stack_data);
%         tmp_is_validQ = ~isoutlier(tmp_fpsd_r);
%         tmp_fpsd_flux = arrayfun(@(x) x.frequency.fpsd_gt_5Hz.('flux').stat.median, stack_data);
        fig_hdl = figure;        
        ax_hdl = axes(fig_hdl);
        bp_hdl = boxplot(ax_hdl, tmp_fpsd_r, is_pia_vessel_Q);
%         ax_hdl.YLim = [1e-4, 1];
        ax_hdl.YLabel.String = sprintf('Fraction of %s power in [%d, %d] Hz', tmp_fn, pps_str.frac_psd_f_range_Hz);
        ax_hdl.YScale = 'log';
        ax_hdl.XTickLabel = {'Penetrating', 'Pia'};
        grid(ax_hdl, 'on');
        
        pen_data = tmp_fpsd_r(is_PA_vessel_Q);
        pia_data = tmp_fpsd_r(is_pia_vessel_Q);
        pen_stat = fun_analysis_get_basic_statistics(pen_data);
        pia_stat = fun_analysis_get_basic_statistics(pia_data);
        [diff_h, diff_p] = kstest2(pen_data, pia_data);
        ax_hdl.Title.String = sprintf('KS-test p-value: %.2e', diff_p);

        fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_frac_pwr_in_%d_to_%d_Hz_%s.png', ...
            dataset_name, int_stack_name, pps_str.frac_psd_f_range_Hz, tmp_fn));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
        
        data_table = struct;
        data_table.fpsd = tmp_fpsd_r;
        data_table.avg_radius_um = vessel_avg_r_list;
        data_table.is_pia = is_pia_vessel_Q;
        data_table = struct2table(data_table);
        table_fp = strrep(fig_fp, 'png', 'csv');
        writetable(data_table, table_fp);
        
%         fig_hdl = figure;
%         ax_hdl = axes(fig_hdl);
%         sc_hdl = scatter(ax_hdl, tmp_fpsd_r(is_pia_vessel_Q & tmp_is_validQ), tmp_fpsd_flux(is_pia_vessel_Q & tmp_is_validQ), ...
%             150, vessel_avg_r_list(is_pia_vessel_Q & tmp_is_validQ), 'o', 'filled');
%         hold(ax_hdl, 'on');
%         sc_hdl_1 = scatter(ax_hdl, tmp_fpsd_r(is_PA_vessel_Q & tmp_is_validQ), tmp_fpsd_flux(is_PA_vessel_Q & tmp_is_validQ), ...
%             150, vessel_avg_r_list(is_PA_vessel_Q & tmp_is_validQ), '*');
%         
%         stat_str = struct;
%         stat_str.pia_fr = mean(tmp_fpsd_r(is_pia_vessel_Q & tmp_is_validQ));
%         stat_str.pia_fr_std = std(tmp_fpsd_r(is_pia_vessel_Q & tmp_is_validQ));
%         stat_str.pia_ff = mean(tmp_fpsd_flux(is_pia_vessel_Q & tmp_is_validQ));
%         stat_str.pia_ff_std = std(tmp_fpsd_flux(is_pia_vessel_Q & tmp_is_validQ));
%         
%         stat_str.pa_ff_a = mean(tmp_fpsd_flux(is_PA_vessel_Q & tmp_is_validQ));
%         stat_str.pa_ff_std = std(tmp_fpsd_flux(is_PA_vessel_Q & tmp_is_validQ));
%         stat_str.pa_ff_a = mean(tmp_fpsd_flux(is_PA_vessel_Q & tmp_is_validQ));
%         stat_str.pa_ff_std = std(tmp_fpsd_flux(is_PA_vessel_Q & tmp_is_validQ));
%         
%         legend_str = {sprintf('Pia\n\t%.3f \\pm %.3f', stat_str.pia_fr, stat_str.pia_fr_std), ...
%             sprintf('PA\n\t%.3f \\pm %.3f', stat_str.pa_ff_a, stat_str.pa_ff_std)};
%         legend(ax_hdl, [sc_hdl, sc_hdl_1], legend_str, 'Location', 'best');
%         
%         ax_hdl.XLim = [0, 0.03];
%         grid(ax_hdl, 'on');
%         ax_hdl.XLabel.String = 'Fractional power for radius';
%         ax_hdl.YLabel.String = 'Fractional power for flux';
%         ax_hdl.YLim = [0,1];
%         sc_hdl.MarkerFaceAlpha = 0.6;
%         sc_hdl_1.MarkerFaceAlpha = 0.6;
%         c_bar = colorbar(ax_hdl);
%         cm_hdl = colormap(ax_hdl, 'jet');
%         c_bar.Label.String = 'Radius (\mum)';
%         ax_hdl.FontSize = 12;
%         ax_hdl.Title.String = sprintf('Fractional power at > 5 Hz');
    end
    %% Radius modulation in frequency band
    tmp_dr2r = arrayfun(@(x) x.radius.r_or.rm_in_fband.stat_mean.mean, stack_data);
%     tmp_is_validQ = ~isoutlier(tmp_dr2r);
    tmp_is_validQ = true(numel(tmp_dr2r), 1);
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
    sc_hdl = scatter(ax_hdl, vessel_avg_r_list(is_pia_vessel_Q & tmp_is_validQ), tmp_dr2r(is_pia_vessel_Q & tmp_is_validQ), ...
        150, 'o', 'filled');
    hold(ax_hdl, 'on');
    sc_hdl_1 = scatter(ax_hdl, vessel_avg_r_list(is_PA_vessel_Q & tmp_is_validQ), tmp_dr2r(is_PA_vessel_Q & tmp_is_validQ), ...
        150, '*');
    ax_hdl.XLabel.String = 'Radius, \langle r \rangle (\mum)';
    ax_hdl.YLabel.String = '\langle|\delta r_H|\rangle / \langle r \rangle ';
    sc_hdl.MarkerFaceAlpha = 0.6;
    sc_hdl_1.MarkerFaceAlpha = 0.6;
    
    pen_data = tmp_dr2r(is_PA_vessel_Q & tmp_is_validQ);
    pia_data = tmp_dr2r(is_pia_vessel_Q & tmp_is_validQ);
    stat_str = struct;
    stat_str.pia_fr = mean(pia_data);
    stat_str.pia_fr_std = std(pia_data);
    stat_str.pa_fr = mean(pen_data);
    stat_str.pa_fr_std = std(pen_data);
    
    [diff_h, diff_p] = kstest2(pen_data, pia_data);
    ax_hdl.Title.String = sprintf('Radius modulation in [%d, %d] Hz\nKS-test p-value: %.2e',  ...
        pps_str.frac_psd_f_range_Hz, diff_p);
    
    legend_str = {sprintf('Pia\n\t%.1e \\pm %.1e', stat_str.pia_fr, stat_str.pia_fr_std), ...
        sprintf('PA\n\t%.1e \\pm %.1e', stat_str.pa_fr, stat_str.pa_fr_std)};
    legend(ax_hdl, [sc_hdl, sc_hdl_1], legend_str, 'Location', 'best');
    
    fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_r_modulation_in_%d_to_%d_Hz_%s.png', ...
        dataset_name, int_stack_name, pps_str.frac_psd_f_range_Hz, tmp_fn));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    
    data_table = struct;
    data_table.fpsd = tmp_dr2r;
    data_table.avg_radius_um = vessel_avg_r_list;
    data_table.is_pia = is_pia_vessel_Q;
    data_table = struct2table(data_table);
    table_fp = strrep(fig_fp, 'png', 'csv');
    writetable(data_table, table_fp);
%% #######################################################################
%% How to summarize the properties of a single ROI? 
plt_group = cell(2, 0);
plt_group(:, end+1) = {'all', true(size(stack_data))};
plt_group(:, end+1) = {'pia', is_pia_vessel_Q};
plt_group(:, end+1) = {'penetrating', is_PA_vessel_Q};
num_plt_group = size(plt_group, 2);
max_vaso_period = 15;

for i_group = 1 : num_plt_group
    tmp_name = plt_group{1, i_group};
    tmp_mask = plt_group{2, i_group};
    flux_vs_r_int_stat = arrayfun(@(x) x.flux.(flux_name).vm_int_stat, ...
        stack_data(tmp_mask), 'UniformOutput', false);
    flux_vs_r_int_stat = cat(1, flux_vs_r_int_stat{:});
    flux_vs_r_int_stat = stack_data.select_vasomition_interval_stat(flux_vs_r_int_stat, ...
        stack_data(1).parameter.vasomotion_min_fraction);
    is_valid_Q = flux_vs_r_int_stat.durations < max_vaso_period;
    flux_vs_r_int_stat = flux_vs_r_int_stat(is_valid_Q, :);
    flux_vs_r_int_stat.dr2r = flux_vs_r_int_stat.r_range ./ flux_vs_r_int_stat.r_mean;
    % flux_vs_r_int_stat = VRRFAI.select_vasomition_interval_stat(flux_vs_r_int_stat, ...
    %     0.33);
    %% R2 vs dr/r
    [fig_hdl, ax_hdl] = VRRFAI.vis_xy_joint_dist_in_hist(flux_vs_r_int_stat.dr2r, ...
        flux_vs_r_int_stat.r2f_r2, 0:0.05:0.7, 0:0.1:1, [], ...
        '\Deltar/\langler\rangle', 'R^2', tmp_name);
    cbar_hdl = colorbar(ax_hdl);
    cbar_hdl.Label.String = 'Number of intervals';
    cbar_hdl.Limits(1) = 0;
    fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_%s_flux_vs_radius_fit_r2_vs_dr2r.png', ...
        dataset_name, int_stack_name, tmp_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    %% k vs dr/r
    [fig_hdl, ax_hdl] = VRRFAI.vis_xy_joint_dist_in_hist(flux_vs_r_int_stat.dr2r, ...
        flux_vs_r_int_stat.dfdr_n_k, 0:0.05:0.7, -2.5:0.25:5, [], ...
        '\Deltar/\langler\rangle', 'dq_n/dr_n', tmp_name);
    cbar_hdl = colorbar(ax_hdl);
    cbar_hdl.Label.String = 'Number of intervals';
    cbar_hdl.Limits(1) = 0;
    fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_%s_flux_vs_radius_fit_k_vs_dr2r.png', ...
        dataset_name, int_stack_name, tmp_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
end
%%
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% histogram2(ax_hdl, flux_vs_r_int_stat.dfdr_n_k, flux_vs_r_int_stat.r2f_r2, ...
%     -4:0.25:4, 0:0.1:1, 'DisplayStyle', 'tile');

fig_hdl = VRRFAI.vis_vasomotion_interval_stat_v2(flux_vs_r_int_stat);
ax_hdl = fig_hdl.Children(2);
ax_hdl.XLim = [-0.5, 3];
sc_hdl = fig_hdl.Children(2).Children;
sc_hdl.MarkerFaceAlpha = 0.2;
sc_hdl.AlphaData = 0.2 * ones(size(sc_hdl.XData));

fig_hdl = VRRFAI.vis_vasomotion_interval_stat_hist2(flux_vs_r_int_stat);
%% Integration analysis
% selected_idx = cell(4, 1);
% selected_idx{1} = [1,2,3,4,5,7,9];
% selected_idx{2} = [1,2,3,4,6,7,11,12,13,14,15,16,17,18];
% selected_idx{3} = [2,3,5];
% selected_idx{4} = [1,2,3,5,6,8,11,12,13,14,15,16,17,19,20,22,23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42];
% 
% merged_data = cellfun(@(x, y) x(y), stack_data, selected_idx, 'UniformOutput', false);
% merged_data = cat(1, merged_data{:});
% merged_data = cat(1, stack_data{:});
% vis_name_prefix = sprintf('%s_TBRLBB', dataset_name);

% stack_idx = 1;
raw_data_cell = arrayfun(@(x) x.raw_data, stack_data, 'UniformOutput', false);
merged_data = cat(1, raw_data_cell{:});
vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    int_stack_name), 'Integrated_single_stack');
vis_name_prefix = int_stack_name;

is_valid_data = arrayfun(@(x) isa(x, 'VRRFAI'), merged_data);
merged_data = merged_data(is_valid_data);
data_label = arrayfun(@(x) sprintf('%sROI%03d%05d', x.stack(7:end), x.roi_idx, x.record_idx),...
    merged_data, 'UniformOutput', false);
    %% Average radius vs average flow     
    r_avg = arrayfun(@(x) x.radius.(r_name).stat.mean, stack_data);
    r_std = arrayfun(@(x) x.radius.(r_name).stat.std, stack_data);
    r_ptl = arrayfun(@(x) x.radius.(r_name).stat.prctile_val(:, [7,8,9]),...
        stack_data, 'UniformOutput', false);
    r_ptl = cat(1, r_ptl{:});
    r_eb = abs(r_ptl(:, [1,3]) - r_ptl(:, 2));
    
    flux_avg = arrayfun(@(x) x.flux.(flux_name).stat.mean, stack_data);
    flux_std = arrayfun(@(x) x.flux.(flux_name).stat.std, stack_data);
    flux_ptl = arrayfun(@(x) x.flux.(flux_name).stat.prctile_val(:, [7,8,9]),...
        stack_data, 'UniformOutput', false);
    flux_ptl = cat(1, flux_ptl{:});
    flux_eb = abs(flux_ptl(:, [1,3]) - flux_ptl(:, 2));
    
    fig_hdl = figure;
    fig_hdl.Position(3) = fig_hdl.Position(3) * 2;
    ax_hdl_1 = subplot(1,2,1);
    for i = 1 : 2
        if i == 1
            errorbar(ax_hdl_1, r_ptl(is_pia_vessel_Q, 2), flux_ptl(is_pia_vessel_Q, 2), ...
                flux_eb(is_pia_vessel_Q, 1), flux_eb(is_pia_vessel_Q, 2), ...
                r_eb(is_pia_vessel_Q, 1), r_eb(is_pia_vessel_Q, 2), 'LineStyle', 'none');
        else
            hold(ax_hdl_1, 'on');
            errorbar(ax_hdl_1, r_ptl(is_PA_vessel_Q, 2), flux_ptl(is_PA_vessel_Q, 2),...
                flux_eb(is_PA_vessel_Q, 1), flux_eb(is_PA_vessel_Q, 2), ...
                r_eb(is_PA_vessel_Q, 1), r_eb(is_PA_vessel_Q, 2), 'LineStyle', 'none');
        end
    end    
    ax_hdl_1.XLabel.String = 'Radius (\mum)';
    ax_hdl_1.YLabel.String = 'Flux (RBC/s)';
    ax_hdl_1.XLim(1) = 0;
    ax_hdl_1.YLim(1) = 0;
    grid(ax_hdl_1, 'on');
    legend(ax_hdl_1, sprintf('#Pia: %d', nnz(is_pia_vessel_Q)),...
        sprintf('#Penetrating: %d', nnz(is_PA_vessel_Q)), ...
        'Location', 'northwest');
    ax_hdl_2 = subplot(1,2,2);
    for i = 1 : 2
        if i == 1
            errorbar(ax_hdl_2, r_ptl(is_pia_vessel_Q, 2), flux_ptl(is_pia_vessel_Q, 2), ...
                flux_eb(is_pia_vessel_Q, 1), flux_eb(is_pia_vessel_Q, 2), ...
                r_eb(is_pia_vessel_Q, 1), r_eb(is_pia_vessel_Q, 2), 'LineStyle', 'none');
        else
            hold(ax_hdl_2, 'on');
            errorbar(ax_hdl_2, r_ptl(is_PA_vessel_Q, 2), flux_ptl(is_PA_vessel_Q, 2),...
                flux_eb(is_PA_vessel_Q, 1), flux_eb(is_PA_vessel_Q, 2), ...
                r_eb(is_PA_vessel_Q, 1), r_eb(is_PA_vessel_Q, 2), 'LineStyle', 'none');
        end
    end    
    ax_hdl_2.XLabel.String = 'Radius (\mum)';
    ax_hdl_2.YLabel.String = 'Flux (RBC/s)';
    ax_hdl_2.XScale = 'log';
    ax_hdl_2.YScale = 'log';
    grid(ax_hdl_2, 'on');
    [linear_fit_hdl] = fun_vis_linear_fit_data_in_loglog(r_ptl(:, 2), ...
        flux_ptl(:, 2), ax_hdl_2);
    
    fig_fp = fullfile(int_fig_folder, sprintf('%s_%s_roi_overall_flux_vs_r.png', ...
        dataset_name, int_stack_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    %% Slow oscillation coherence peak 
    coh_group_name = {'Radius', 'Flux'};
    num_rois = numel(merged_data);
    test_data = merged_data(1);
    slow_frequency_range = [0.05, 0.25];
    frequency_list = test_data.frequency.r.psd.f;
    is_slow_frequency_Q = frequency_list >= slow_frequency_range(1) & ...
        frequency_list <= slow_frequency_range(2);
    low_frequency_list = frequency_list(is_slow_frequency_Q);
    %%
    coh_pair_idx_1 = 1;
    coh_pair_idx_2 = 2;
    
    coh_pair_name_1 = coh_group_name{coh_pair_idx_1};
    coh_pair_name_2 = coh_group_name{coh_pair_idx_2};
    
    [coh_peak_f, coh_peak, coh_peak_phase] = deal(nan(num_rois, 1));
    for iter_file = 1 : num_rois
        tmp_data = merged_data(iter_file);
        tmp_coh_abs = tmp_data.frequency.coherence.abs{coh_pair_idx_1, coh_pair_idx_2}(is_slow_frequency_Q);
        tmp_coh_phase = tmp_data.frequency.coherence.arg{coh_pair_idx_1, coh_pair_idx_2}(is_slow_frequency_Q);
        [tmp_coh_peak, tmp_peak_idx] = max(tmp_coh_abs);
        coh_peak_f(iter_file) = low_frequency_list(tmp_peak_idx);
        coh_peak(iter_file) = tmp_coh_peak;
        coh_peak_phase(iter_file) = tmp_coh_phase(tmp_peak_idx);
    end
    coh_peak_phase_d = coh_peak_phase / pi * 180;
    
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
    ax_hdl_1 = subplot(2,2,1);
    scatter(ax_hdl_1, coh_peak_f, coh_peak);
    ax_hdl_1.XLabel.String = 'Frequency (Hz)';
    ax_hdl_1.YLabel.String = 'Peak coherence';
    ax_hdl_1.YLim = [0,1];
    grid(ax_hdl_1, 'on');
    
    ax_hdl_2 = subplot(2,2,2);
    scatter(ax_hdl_2, coh_peak_f, coh_peak_phase_d);
    ax_hdl_2.XLabel.String = 'Frequency (Hz)';
    ax_hdl_2.YLabel.String = 'Peak phase (^\circ)';
    ax_hdl_2.YLim = [-180, 180];
    grid(ax_hdl_2, 'on');
    
    ax_hdl_3 = subplot(2,2,3);
    scatter(ax_hdl_3, coh_peak, coh_peak_phase_d);
    ax_hdl_3.XLabel.String = 'Peak coherence';
    ax_hdl_3.XLim = [0, 1];
    ax_hdl_3.YLabel.String = 'Peak phase (^\circ)';
    ax_hdl_3.YLim = [-180, 180];
    grid(ax_hdl_3, 'on');
    
    ax_hdl_5 = subplot(2,2,4);
    p_d_stat = fun_analysis_get_basic_statistics(coh_peak_phase_d);
    [test_h, test_p] = ttest(coh_peak_phase_d);
    histogram(ax_hdl_5, coh_peak_phase_d, 'Normalization', 'count');
    ax_hdl_5.XLabel.String = 'Peak Phase (^\circ)';
    ax_hdl_5.XLim = [-180, 180];
    ax_hdl_5.YLabel.String = 'Count';
    grid(ax_hdl_5, 'on');
    legend(ax_hdl_5, sprintf('%.1f \\pm %.1f\nT-test p: %.2e', ...
        p_d_stat.mean, p_d_stat.std, test_p));
    
    fig_fp = fullfile(vis_folder, sprintf('%s_%s', vis_name_prefix,...
        sprintf('%s_vs_%s_coherence_peak_scatters.png', coh_pair_name_1, coh_pair_name_2)));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    
    %% 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
