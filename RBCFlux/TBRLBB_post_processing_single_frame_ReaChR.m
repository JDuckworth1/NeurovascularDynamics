clc;clear;close all;
DataManager = FileManager;
dataset_name = 'DKLab';
%stack_list = {'TBRLBB20211129', 'TBRLBB20211119', 'TBRLBB20211118',...
    %'TBRLBB20211103', 'TBRLBB20211108'};
stack_list = {'JD230523R1495trials_011524'}; %'20231219R155E'
num_stack = numel(stack_list);
dev_mode = false;
% Processing parameters
ppp_str = struct;
ppp_str.outlier_multiplifier = 3;
ppp_str.outlier_window_s = 1;
ppp_str.sm_wd_s = 2;
% Peak detection
ppp_str.min_dilation_fraction = 0.03;
ppp_str.min_dilation_spacing_s = 2;
% Binning
ppp_str.bin_num_phase = 16;
ppp_str.bin_num_radius = 16;
ppp_str.bin_num_rbc_count = 16;

ppp_str.feq.f_resolve_Hz = 0.05;
ppp_str.feq.fft_pad_factor = 2;
ppp_str.feq.vasomotion_feq_Hz = [0.05, 0.25];


ppp_str.speed.rbc_l_um = 5;
ppp_str.speed.t_std_coeff = 3;

vis_opt = struct;
vis_opt.generate_video_Q = false;

%%
for stack_idx = 1 : num_stack
    %%
    stack = stack_list{stack_idx};
    data_info = DataManager.load_data(fullfile(DataManager.fp_raw_data_folder(dataset_name, ...
        stack), sprintf('%s_%s_raw_data_info.mat', dataset_name, stack)));
    
     for iter_file = 1 : data_info.num_files

        roi_idx = data_info.roi_list(iter_file);
        record_idx = data_info.record_idx_list(iter_file);
        %% Load data
        tic_vis = tic;
        fprintf('Processing ROI %d record %d\n', roi_idx, record_idx);
        result_fp = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack, ...
            sprintf('%s_%s_roi_%03d_%05d_sf.mat', dataset_name, stack,...
            roi_idx, record_idx)));
        %result_fp = strrep(result_fp,'analysis_data','analysis_data_fullR');
        result_fp = strrep(result_fp,'analysis_data','analysis_data_97frac');
        
        try
            result_str = DataManager.load_data(result_fp);
            if dev_mode
                tmp_im_fp = fullfile(DataManager.fp_raw_data_folder(dataset_name, stack), ...
                    result_str.radius_data_fn);
                vsl_im = DataManager.load_data(tmp_im_fp);
                vsl_im = vsl_im.im;
            end
        catch ME
            warning('Failed to load data for ROI %d', roi_idx);
            continue;
        end
        
        vis_data = VRRFA(result_str);
        vis_data.parse_parameters(ppp_str);
        try
            %% Process traces
            % Radius
            
%             vis_data.radius.r_or = vis_data.process_radius_data([vis_data.data.radius.radius_um].');
            % Use the radius estimation from each line scan. Reshape the
            % trace into num_line_per_frame x num_frames array, compute the
            % median. 
            radius_l_data = reshape(vis_data.data.radius_l_um, vis_data.parameters.num_line_per_frame, []);
            radius_l_data_med = median(radius_l_data, 1).';
            vis_data.radius.r_or = vis_data.process_radius_data(radius_l_data_med);
            
            % Flux
            rbc_count_in_bin = vis_data.compute_bin_avg_in_windows(...
                vis_data.data.flux.rbc_trace, vis_data.parameters.num_line_per_frame).' * ...
                vis_data.parameters.linescan_rate_Hz;
            vis_data.flux.per_frame = vis_data.process_flux_data(rbc_count_in_bin);          
            % RBC passing time

            rbs_v = vis_data.estimate_rbc_speed(vis_data.parameters.speed.rbc_l_um, ...
                vis_data.parameters.speed.t_std_coeff);            
            vis_data.speed.per_frame = vis_data.process_flux_data(rbs_v);
            %% Compute smoothed signals
                %% long pass filder
%             r_lp = lowpass(vis_data.radius.r_or.data, vis_data.parameters.feq.vasomotion_feq_Hz(2)/2, ...
%                 vis_data.parameters.frame_rate_Hz);
%             f_lp = lowpass(vis_data.flux.per_frame.data, vis_data.parameters.feq.vasomotion_feq_Hz(2)/2, ...
%                 vis_data.parameters.frame_rate_Hz);
%             vis_data.radius.r_or_lp = vis_data.process_radius_data(r_lp);
%             vis_data.flux.per_frame_lp = vis_data.process_flux_data(f_lp);
%             fig_hdl = vis_data.vis_radius_and_flux_with_r_phase(vis_data.parameters.t_s, ...
%                 vis_data.radius.r_or_lp.data, vis_data.flux.per_frame_lp.data,...
%                 vis_data.radius.r_or_lp.phase_l);
                %% Moving mean / median?
            r_sm = movmean(vis_data.radius.r_or.data, vis_data.parameters.sm_wd_sz);
            f_sm = movmean(vis_data.flux.per_frame.data, vis_data.parameters.sm_wd_sz);
            rbc_v_sm = movmean(vis_data.speed.per_frame.data, vis_data.parameters.sm_wd_sz);
            
            vis_data.radius.r_or_mm = vis_data.process_radius_data(r_sm);
            vis_data.flux.per_frame_mm = vis_data.process_flux_data(f_sm);            
            vis_data.speed.per_frame_mm = vis_data.process_flux_data(rbc_v_sm);
            %% Visualize smoothed traces
            % Flux
            fig_hdl = vis_data.vis_radius_vs_speed(vis_data.parameters.t_s, ...
                vis_data.radius.r_or_mm.data, vis_data.flux.per_frame_mm.data, ...
                'Radius (\mum)', 'Flux (/s)', '$dr/\bar{r}$', '$dq/\bar{q}$');
            vis_data.save_figure(fig_hdl, {'average', 'traces'}, ...
                'smoothed_r_f_fit.png', ~dev_mode);    
            % FOR SUPPLEMENT
            r_trace_frame = vis_data.radius.r_or.data;
            r_trace = vis_data.radius.r_or_mm.data;
            t_trace = vis_data.parameters.t_s;
            flux_trace = vis_data.flux.per_frame_mm.data;
            pctlim = 0.20;
            r_20pct = r_trace(1) - pctlim*r_trace(1);
            r_diff = r_trace - r_20pct;
            rmin_vals = find(r_diff < 0); %Use first passing
            r20pt = rmin_vals(1);
            rtrace_20pt = r_trace_frame(1:r20pt);
            ttrace_20pt = t_trace(1:r20pt);
            fluxtrace_20pt = flux_trace(1:r20pt);
            figure; plot(ttrace_20pt,2*rtrace_20pt); hold on; yyaxis right;
            plot(ttrace_20pt,fluxtrace_20pt);
            ax = gca; ax.TickLabelInterpreter = 'latex';
            print(gcf, '-depsc2',['JD230523R1495trials_011524_Line10_Trial1_20pt_Traces'])
            figure; plot(t_trace,2*r_trace_frame); hold on; yyaxis right;
            plot(t_trace,flux_trace);
            ax = gca; ax.TickLabelInterpreter = 'latex';
            print(gcf, '-depsc2',['JD230523R1495trials_011524_Line10_Trial1_full_Traces'])
            
            % Speed
            fig_hdl = vis_data.vis_radius_vs_speed(vis_data.parameters.t_s, ...
                vis_data.radius.r_or_mm.data, vis_data.speed.per_frame_mm.data, ...
                'Radius (\mum)', 'Speed (mm/s)', '$dr/\bar{r}$', '$dv/\bar{v}$');
            vis_data.save_figure(fig_hdl, {'average', 'traces'}, ...
                'smoothed_r_v_fit.png', ~dev_mode);
            %% Vasodilation intervals
            [vis_data.radius.r_or_mm.vm_int, vis_data.radius.r_or_mm.phase_l, fig_hdl] = ...
                vis_data.compute_vasodilation_intervals(...
                vis_data.radius.r_or_mm.data, vis_data.parameters.min_dilation_fraction, ...
                vis_data.parameters.min_dilation_spacing_s, true);
            
            vis_data.save_figure(fig_hdl, {'average', 'traces'},...
                'r_or_mm_w_phase_l.png', ~dev_mode);
                   %% Visualize radius oscillation phase
            fig_hdl = vis_data.vis_radius_and_flux_with_r_phase(vis_data.parameters.t_s, ...
                vis_data.radius.r_or_mm.data, vis_data.flux.per_frame_mm.data,...
                vis_data.radius.r_or_mm.phase_l);        
            vis_data.save_figure(fig_hdl, {'average', 'traces'}, ...
                'smoothed_r_f_rphase_l.png', ~dev_mode);
            %% Analyze intervals
            [vis_data.flux.per_frame_mm.vm_int, vis_data.flux.per_frame_mm.vm_int_stat] =...
                vis_data.analyze_vasomotion_interval(vis_data.radius.r_or_mm.vm_int, ...
                vis_data.radius.r_or_mm.data, vis_data.flux.per_frame_mm.data);
            
            [vis_data.speed.per_frame_mm.vm_int, vis_data.speed.per_frame_mm.vm_int_stat] = ...
                vis_data.analyze_vasomotion_interval(vis_data.radius.r_or_mm.vm_int, ...
                vis_data.radius.r_or_mm.data, vis_data.speed.per_frame_mm.data);
            %% Visualize flux - radius phase scatter plot
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            histogram2(ax_hdl, vis_data.radius.r_or_mm.phase_l, vis_data.flux.per_frame_mm.data, ...
                16, ...
                'DisplayStyle', 'tile');
            ax_hdl.YLabel.String = 'Flux (RBC/s)';
            ax_hdl.XLabel.String = 'Radius phase (rad)';
            
            vis_data.save_figure(fig_hdl, {'vs_phase'}, 'f_vs_pl_all.png', ...
                ~dev_mode);               
            %% Flux vs estiamted flux
            vis_data.speed.per_frame_mm.est_flux = vis_data.speed.per_frame_mm.data .* ...
                (vis_data.radius.r_or_mm.data .^ 2) * pi / 1e6;
            vis_data.speed.per_frame_mm.est_flux_stat = fun_analysis_get_basic_statistics(...
                vis_data.speed.per_frame_mm.est_flux);
            vis_data.speed.per_frame_mm.est_flux_n = (vis_data.speed.per_frame_mm.est_flux - ...
                vis_data.speed.per_frame_mm.est_flux_stat.mean) ./ ...
                vis_data.speed.per_frame_mm.est_flux_stat.mean;
            
            fig_hdl = vis_data.vis_radius_vs_speed(vis_data.parameters.t_s, ...
                vis_data.speed.per_frame_mm.est_flux, ...
                vis_data.flux.per_frame_mm.data, ...                
                'Volume Flux Q (nL/s)', 'Flux q (RBC/s)', ...
                '$dQ/\bar{Q}$', '$dq/\bar{q}$');
            vis_data.save_figure(fig_hdl, {'flux_rbc_vs_volume'}, ...
                'f_vs_vr2_all.png', ~dev_mode);
            %% Spectrum analysis
            vis_data.frequency = vis_data.compute_frequency_analysis(...
                vis_data.radius.r_or.data, vis_data.flux.per_frame.data, ...
                vis_data.speed.per_frame.data);
            % Plot power spectrum density
            fig_hdl = vis_data.vis_power_spectrum_densities(...
                {vis_data.frequency.r.psd, vis_data.frequency.flux.psd, vis_data.frequency.speed.psd}, ...
                {'Radius', 'Flux', 'Speed'});
            vis_data.save_figure(fig_hdl, {'average', 'PSD'}, ...
                'radius_flux_speed_SPD.png', ~dev_mode);
            vis_data.vis_pairwise_coherence(~dev_mode);
            %%
            vis_data.data = [];
            DataManager.write_data(vis_data.filepath, vis_data);

            
            fprintf('Finish visualizing data for ROI %d. Elapsed time is %.2f seconds\n', ...
                roi_idx, toc(tic_vis));
        catch ME
            fprintf('Fail to visualize %s. Skip.\n', vis_data.vis_fp_template);
            rethrow(ME);
            %             fprintf('Error message: %s\n', ME.message);
        end
    end
end