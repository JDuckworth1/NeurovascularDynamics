classdef VRRFA < handle
    % Vessel radius and RBC flux - analysis
    properties
        data struct
        dataset_name char
        stack char
        roi_idx (1,1) double
        record_idx (1,1) double
        data_fp char
        vis_fp_template char
        %         vis_root_folder char
        %         filepath char
        
        parameters
        % Results
        radius
        radius_l
        flux
        speed
        frequency
        %
        
    end
    
    
    
    methods
        function obj = VRRFA(result_str)
            
            % Intialization
            obj.data = result_str;
            obj.dataset_name = result_str.dataset_name;
            obj.stack = result_str.stack;
            obj.roi_idx = result_str.roi_idx;
            obj.record_idx = result_str.record_idx;
            obj.vis_fp_template = sprintf('%s_%s_ROI_%03d_%05d', obj.dataset_name, ...
                obj.stack, obj.roi_idx, obj.record_idx);
        end
        
        function fp = vis_root_folder(obj)
            DataManager = FileManager;
            fp = sprintf('%s_sl', DataManager.fp_visualization_folder(obj.dataset_name, ...
                obj.stack));
        end
        
        function fp = filepath(obj)
            DataManager = FileManager;
            fp = fullfile(DataManager.fp_analysis_data_folder(obj.dataset_name, ...
                obj.stack), 'vis_data', sprintf('%s_vis_data.mat', obj.vis_fp_template));
        end
        function fp = filepath_95(obj)
            DataManager = FileManager;
            fp = fullfile(DataManager.fp_analysis_data_folder(obj.dataset_name, ...
                obj.stack), 'vis_data', sprintf('%s_vis_data.mat', obj.vis_fp_template));
            fp = strrep(fp,'analysis_data','analysis_data_3_4R_95frac');
        end        
        function fp = filepath_97(obj)
            DataManager = FileManager;
            fp = fullfile(DataManager.fp_analysis_data_folder(obj.dataset_name, ...
                obj.stack), 'vis_data', sprintf('%s_vis_data.mat', obj.vis_fp_template));
            fp = strrep(fp,'analysis_data','analysis_data_3_4R_97frac');
        end
        function fp = filepath_98(obj)
            DataManager = FileManager;
            fp = fullfile(DataManager.fp_analysis_data_folder(obj.dataset_name, ...
                obj.stack), 'vis_data', sprintf('%s_vis_data.mat', obj.vis_fp_template));
            fp = strrep(fp,'analysis_data','analysis_data_3_4R_98frac');
        end
        function fp = filepath_99(obj)
            DataManager = FileManager;
            fp = fullfile(DataManager.fp_analysis_data_folder(obj.dataset_name, ...
                obj.stack), 'vis_data', sprintf('%s_vis_data.mat', obj.vis_fp_template));
            fp = strrep(fp,'analysis_data','analysis_data_3_4R_99frac');
        end        
        
        function save_figure(obj, fig_hdl, folder_cell, fn_postfix, delete_hdl_Q)
            if nargin < 5
                delete_hdl_Q = false;
            end
            
            if iscell(folder_cell)
                vis_folder = fullfile(obj.vis_root_folder, folder_cell{:});
            else
                vis_folder = fullfile(obj.vis_root_folder, folder_cell);
            end
            fig_fp = fullfile(vis_folder, sprintf('%s_%s', obj.vis_fp_template, ...
                fn_postfix));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            if delete_hdl_Q
                delete(fig_hdl);
            end
        end
        
        function obj = parse_parameters(obj, ppp_str)
            obj.parameters = ppp_str;
            
            obj.parameters.num_frame = obj.data.parameters.num_frame;
            obj.parameters.frame_rate_Hz = obj.data.parameters.frame_rate_Hz;
            obj.parameters.dt_s = obj.data.parameters.dt_s;
            obj.parameters.total_t_s = obj.parameters.num_frame * obj.parameters.dt_s;
            obj.parameters.linescan_rate_Hz = obj.data.parameters.linescan_rate_Hz;
            obj.parameters.num_line_per_frame = obj.data.parameters.im_size(2);
            
            obj.parameters.outlier_window = round(ppp_str.outlier_window_s / 2 * ...
                obj.parameters.frame_rate_Hz) * 2 + 1;
            obj.parameters.sm_wd_sz = round(ppp_str.sm_wd_s / 2 * ...
                obj.parameters.frame_rate_Hz) * 2 + 1;
            
            obj.parameters.feq.f_Rayleigh_Hz = 1 ./ obj.parameters.total_t_s;
            obj.parameters.feq.p = ceil(obj.parameters.feq.f_resolve_Hz / ...
                obj.parameters.feq.f_Rayleigh_Hz / 2);
            obj.parameters.t_s = (1 : obj.parameters.num_frame) * obj.parameters.dt_s;
            
        end
        
        function flux_str = process_flux_data(obj, rbc_count_trace)
            flux_str = struct;
            flux_str.data = rbc_count_trace;
            % Fill missing value
            flux_str.data = fillmissing(flux_str.data, 'linear');
            flux_str.stat = fun_analysis_get_basic_statistics(flux_str.data);
        end
        
        function rbc_v_mm_s = estimate_rbc_speed(obj, est_rbc_l_um, ...
                passing_std_coeff)
            line_idx_edge = 0 : obj.data.parameters.num_lines_per_frame : ...
                obj.data.parameters.num_lines;
            rbc_idx = fun_bin_data_to_idx_list_by_edges(obj.data.flux.t_mean, ...
                line_idx_edge);
            rbc_t_std = fun_bin_data_to_cells_by_ind(obj.data.flux.t_std, ...
                rbc_idx);
            rbc_t_std_mean = cellfun(@mean, rbc_t_std);            
            
            rbc_v_mm_s = est_rbc_l_um * 1e-3 ./(rbc_t_std_mean * passing_std_coeff / ...
                obj.data.parameters.linescan_rate_Hz);            
        end
        
        function obj = smooth_traces(obj, sm_wd_s)
            % Smooth window size
            sm_wd = round(sm_wd_s * obj.parameters.frame_rate_Hz);
            
            if sm_wd > 1
                obj.radius.radius_um_sm = movmean(obj.radius.radius_um, sm_wd);
                obj.flux.per_frame_sm = movmean(obj.flux.per_frame, sm_wd);
            else
                obj.radius.radius_um_sm = obj.radius.radius_um;
                obj.flux.per_frame_sm = obj.flux.per_frame;
            end
        end
        
        
        function radius_str = process_radius_data(obj, radius_vec)
            radius_str = struct;
            radius_str.data = radius_vec;
            % Reject outliers by percentile - maybe need local stat
            % estimation, otherwise, vessels are actively oscillating
            radius_str.is_outlier_Q = isoutlier(radius_vec, ...
                'movmean', obj.parameters.outlier_window, 'ThresholdFactor', ...
                obj.parameters.outlier_multiplifier);
            radius_str.frac_outlier = mean(radius_str.is_outlier_Q);
            if any(radius_str.is_outlier_Q)
                fprintf('Fraction of radius outliers: %.2e. Replace by linear interpolation\n', ...
                    radius_str.frac_outlier);
            else
                fprintf('All radius measurements are inliers\n');
            end
            radius_str.data(radius_str.is_outlier_Q) = nan;
            radius_str.data = fillmissing(radius_str.data, 'linear');
            
            radius_str.stat = fun_analysis_get_basic_statistics(radius_str.data);
            radius_str.phase = obj.compute_r_trace_phase_ht(radius_str.data);
        end        
        
        function [result_str, stat_str] = analyze_vasomotion_interval(obj, vd_interval, ...
                r_trace, y_trace)
            num_interval = numel(vd_interval);
            
            result_str = cell(num_interval, 1);            
            for i = 1 : num_interval
                tmp_itv = vd_interval(i);
                tmp_ind = tmp_itv.x_ind;
                tmp_itv.r = r_trace(tmp_ind);                
                tmp_itv.y = y_trace(tmp_ind);                           
                
                tmp_itv.r_stat = fun_analysis_get_basic_statistics(tmp_itv.r);
                tmp_itv.y_stat = fun_analysis_get_basic_statistics(tmp_itv.y);
                
                tmp_itv.r_n = (tmp_itv.r - tmp_itv.r_stat.mean) / tmp_itv.r_stat.mean;
                tmp_itv.y_n = (tmp_itv.y - tmp_itv.y_stat.mean) / tmp_itv.y_stat.mean;
                
                tmp_itv.r2f_lc = obj.compute_linear_regression_parameter(tmp_itv.r,...
                    tmp_itv.y);
                tmp_itv.dfdr_n_k = tmp_itv.r2f_lc(1) *...
                    tmp_itv.r_stat.mean / tmp_itv.y_stat.mean;   
                tmp_corr = corrcoef(tmp_itv.r, tmp_itv.y);
                tmp_itv.r2f_r2 = tmp_corr(2) ^ 2;
                % Fit with cos?                
                result_str{i} = tmp_itv;
            end
            result_str = cat(1, result_str{:});
            
            stat_str = struct;
            stat_str.durations = [result_str.duration_s].';
            stat_str.rf_frac = cat(1, result_str.rf_frac);
            stat_str.r_mean = arrayfun(@(x) x.r_stat.mean, result_str);
            stat_str.r_min = arrayfun(@(x) x.r_stat.min, result_str);
            stat_str.r_max = arrayfun(@(x) x.r_stat.max, result_str);
            stat_str.r_range = stat_str.r_max - stat_str.r_min;
            stat_str.y_mean = arrayfun(@(x) x.y_stat.mean, result_str);
            stat_str.dfdr_n_k = [result_str.dfdr_n_k].';
            stat_str.r2f_r2 = [result_str.r2f_r2].';
            stat_str.rf_s = stat_str.rf_frac .* stat_str.durations;
            stat_str = struct2table(stat_str);
        end
        
        function feq_str = compute_frequency_analysis(obj, r_trace, flux_trace, ...
                speed_trace)
            feq_str = struct;
            feq_str.r = obj.compute_power_spectrum_density(r_trace, ...
                obj.parameters.feq.p, obj.parameters.feq.fft_pad_factor, obj.parameters.frame_rate_Hz);
            feq_str.flux = obj.compute_power_spectrum_density(flux_trace, ...
                obj.parameters.feq.p, obj.parameters.feq.fft_pad_factor, obj.parameters.frame_rate_Hz);
            feq_str.speed = obj.compute_power_spectrum_density(speed_trace, ...
                obj.parameters.feq.p, obj.parameters.feq.fft_pad_factor, obj.parameters.frame_rate_Hz);
            
            feq_str.coherence = fun_analysis_feq_compute_coherence_dpss(...
                cat(3, feq_str.r.taper_estimates, feq_str.flux.taper_estimates, feq_str.speed.taper_estimates));
        end
        
        function obj = analyze_slow_cofluctuation(obj, r_trace, flux_trace, mov_wd_sz, save_fig_Q)
            if nargin < 6
                save_fig_Q = false;
            end
            vs_flow = struct;
            vs_flow.r_mm = movmean(r_trace, mov_wd_sz, 'omitnan');
            vs_flow.flux_mm = movmean(flux_trace, mov_wd_sz, 'omitnan');
            %% Compute traces
            field_name = fieldnames(vs_flow);
            for iter_fn = 1 : numel(field_name)
                tmp_fn = field_name{iter_fn};
                tmp_data = vs_flow.(tmp_fn);
                tmp_stat = fun_analysis_get_basic_statistics(tmp_data);
                tmp_df2f = (tmp_data - tmp_stat.mean) ./ tmp_stat.mean;
                vs_flow.(sprintf('%s_stat', tmp_fn)) = tmp_stat;
                vs_flow.(sprintf('df2f_%s', tmp_fn)) = tmp_df2f;
            end
            obj.radius.vs_flow = vs_flow;
            %% Compute cofluctuation
            trace_r = obj.radius.vs_flow.df2f_r_mm;
            x_name = 'Radius';
            plot_cell = cell(2, 0);
            plot_cell(:, end+1) = {obj.radius.vs_flow.df2f_flow_mm, 'Flow'};
            plot_cell(:, end+1) = {obj.radius.vs_flow.df2f_flux_mm, 'Flux'};
            
            num_pairs = size(plot_cell, 2);
            for iter_pair = 1 : num_pairs
                tmp_trace_y = plot_cell{1, iter_pair};
                tmp_y_name = plot_cell{2, iter_pair};
                tmp_cov_str = obj.compute_multiscale_moving_trace_correlation(trace_r, ...
                    tmp_trace_y, obj.parameters.cov.wd_list, ...
                    obj.parameters.cov.wd_step_list, obj.parameters.cov.sm_diff_gf_sigma);
                tmp_cov_str.window_s_list = obj.parameters.cov.wd_list_s;
                
                obj.radius.vs_flow.(sprintf('cov_%s', lower(tmp_y_name))) = tmp_cov_str;
                
                fig_hdl = obj.vis_multiscale_moving_correlation(trace_r, tmp_trace_y, tmp_cov_str, ...
                    obj.parameters.dt_s, {x_name, tmp_y_name});
                if save_fig_Q
                    fig_fp = fullfile(obj.vis_root_folder, 'average', 'cov', tmp_y_name, ...
                        sprintf('%s_%s_%s_cov.png', obj.vis_fp_template, x_name, ...
                        tmp_y_name));
                    fun_print_image_in_several_formats(fig_hdl, fig_fp);
                    delete(fig_hdl);
                end
            end
            
        end
        
        function [vd_int_c, phase_lin, varargout] = compute_vasodilation_intervals(...
                obj, r_trace_sm, min_dilation_fraction, min_spacing_s, vis_Q)
            if nargin < 5
                vis_Q = false;
            end
            sample_rate = obj.parameters.frame_rate_Hz;
            avg_r = mean(r_trace_sm, 'omitnan');
            [r_peaks, peak_ind] = findpeaks(r_trace_sm, ...
                'MinPeakProminence', min_dilation_fraction * avg_r, ...
                'MinPeakDistance', sample_rate * min_spacing_s);
            [r_min, min_ind] = findpeaks(-r_trace_sm, ...
                'MinPeakProminence', min_dilation_fraction * avg_r, ...
                'MinPeakDistance', sample_rate * min_spacing_s);
            r_min = -r_min;
            
            if vis_Q
                fig_hdl = figure;
                ax_hdl = axes(fig_hdl);
                plot(ax_hdl, obj.parameters.t_s, r_trace_sm);
                hold(ax_hdl, 'on');
                scatter(ax_hdl, peak_ind / sample_rate, r_peaks, ...
                    'o', 'filled');
                scatter(ax_hdl, min_ind / sample_rate, r_min, ...
                    'o', 'filled');
                ax_hdl.XLabel.String = 'Time (s)';
                ax_hdl.YLabel.String = 'Radius (\mum)';
                yyaxis(ax_hdl, 'right');
                ax_hdl.YLabel.String = 'Phase (rad)';
                ax_hdl.YLim = [-pi, pi];
            end
            % Find vasomotion intervals
            extrema_ind_0 = cat(1, min_ind, peak_ind);
            [extrema_ind, sort_ind] = sort(extrema_ind_0, 'ascend');
            is_peak_Q = sort_ind > numel(min_ind);
            num_peak = numel(peak_ind);
            num_valid_peak = 0;
            vd_int_c = cell(num_peak, 1);
            for i = 1 : (numel(is_peak_Q) - 1)
                if is_peak_Q(i) && i > 1
                    if ~is_peak_Q(i-1) && ~is_peak_Q(i+1)
                        num_valid_peak = num_valid_peak + 1;
                        tmp_ind = extrema_ind((i-1): (i+1)).';
                        tmp_ind(3) = tmp_ind(3) - 1;
                        % Linear approximate of the phase
                        tmp_x_1 = tmp_ind(1) : tmp_ind(2);
                        tmp_p_1 = linspace(-pi, 0, tmp_ind(2) - tmp_ind(1) + 1);
                        tmp_x_2 = (tmp_ind(2) + 1) : tmp_ind(3);
                        tmp_p_2 = linspace(0, pi, tmp_ind(3) - tmp_ind(2) + 1);
                        
                        tmp_str = struct;
                        tmp_str.ep_ind = tmp_ind;
                        tmp_str.ep_ind_l = tmp_ind - tmp_ind(1) + 1;
                        tmp_str.width = tmp_ind(3) - tmp_ind(1);
                        tmp_str.duration_s = tmp_str.width / sample_rate;
                        tmp_str.rf_frac = diff(tmp_ind) / tmp_str.width;
                        
                        tmp_str.x_ind = cat(1, tmp_x_1', tmp_x_2');
                        tmp_str.t_s = tmp_str.x_ind / sample_rate;
                        tmp_str.phase = cat(1, tmp_p_1', tmp_p_2(2:end)');
                        vd_int_c{num_valid_peak} = tmp_str;
                        if vis_Q
                            plot(ax_hdl, tmp_str.t_s, tmp_str.phase, 'g-');
                            line(ax_hdl, tmp_str.t_s([1,1]), ax_hdl.YLim, 'LineStyle', '-.');
                            line(ax_hdl, tmp_str.t_s([end,end]), ax_hdl.YLim, 'LineStyle', '-.');
                        end
                    end
                end
            end
            vd_int_c = vd_int_c(1 : num_valid_peak);
            vd_int_c = cat(1, vd_int_c{:});
            
            phase_lin = nan(size(r_trace_sm));
            phase_lin(cat(1, vd_int_c.x_ind)) = cat(1, vd_int_c.phase);
            
            if nargout > 2
                varargout{1} = fig_hdl;
            end
        end
    end
    %% Process raw traces
    methods(Static)
        function r_phase = compute_r_trace_phase_ht(r_trace, fpass, fs)
            % fpass: bandpass filter frequency range, Hz
            % fs: sampling frequency, scalar, Hz
            if nargin > 2
                assert(numel(fpass) == 2 & isscalar(fs));
                r_trace = bandpass(r_trace, fpass, fs);
            end
            
            r_ht = hilbert(r_trace - mean(r_trace));
            r_phase = angle(r_ht);
        end
        
        function feq_str = compute_power_spectrum_density(data, p, fft_pad_factor, sample_rate_Hz)
            feq_str = struct;
            assert(isvector(data));
            if isrow(data)
                data = data .';
            end
            is_finite_Q = isfinite(data);
            finite_ratio = mean(is_finite_Q);
            if finite_ratio ~= 1
                warning('%.2f%% of the time trace is not finite. Replace by linear interpolation\n', ...
                    (1 - finite_ratio) * 100);
                data(~is_finite_Q) = nan;
                data = fillmissing(data, 'linear');
            end
            [feq_str.psd, feq_str.taper_estimates] = fun_analysis_feq_compute_psd_dpss(...
                data, p, fft_pad_factor, sample_rate_Hz);
            % Compute cumulative power:
            feq_str.psd.cum_pwr = cumsum(feq_str.psd.psd);
            feq_str.psd.cum_pwr_n = feq_str.psd.cum_pwr ./ feq_str.psd.cum_pwr(end);
        end
        
        
        function r_trace_avg = compute_average_r_trace(r_trace, is_valid_x_Q, ...
                r_sampled_pos_x, min_valid_fraction)
            
            % r_trace: num_pos x num_frame
            % is_valid_x_Q: num_col x num_frame
            % r_sampled_pos_x: num_pos x 1
            % min_valid_fraction: scalar in [0, 1]
            
            valid_fraction = mean(is_valid_x_Q, 2);
            valid_fraction = valid_fraction(r_sampled_pos_x);
            is_valid_trace_Q = valid_fraction > min_valid_fraction;
            r_trace = r_trace(is_valid_trace_Q, :);
            r_trace = fun_DKLab_TBRL_fill_missing_radius(r_trace);
            r_trace_avg = mean(r_trace, 1, 'omitnan').';
        end
        
        function avg_vec = compute_bin_avg_in_windows(data_vec, wd_sz)
            num_pts = floor(numel(data_vec) / wd_sz);
            data_vec = data_vec(1:(num_pts * wd_sz));
            data_vec = reshape(data_vec, wd_sz, []);
            avg_vec = mean(data_vec, 1, 'omitnan');
        end
    end
    %% Analyze slow co-fluctuation
    
    methods
        function mt_corr_str = compute_multiscale_moving_trace_correlation(obj, trace_1, trace_2, ...
                wd_sz_list, wd_step_list, gf_list)
            
            if nargin < 6
                gf_list = [];
            end
            
            if isscalar(wd_step_list)
                wd_step_list = repelem(wd_step_list, numel(wd_sz_list), 1);
            end
            
            if isscalar(gf_list)
                gf_list = repelem(gf_list, numel(wd_sz_list), 1);
            end
            
            mt_corr_str = struct;
            mt_corr_str.window_size_list = wd_sz_list;
            mt_corr_str.window_step_list = wd_step_list;
            mt_corr_str.num_windows = numel(wd_sz_list);
            
            num_points = numel(trace_1);
            [rho_mat, rho_p_mat, rho_d_mat, rho_d_p_mat, k_mat, b_mat] = deal(nan(num_points, mt_corr_str.num_windows));
            
            for iter_wd = 1 : mt_corr_str.num_windows
                if ~isempty(gf_list)
                    tmp_corr_str = obj.compute_moving_trace_correlation(...
                        trace_1, trace_2, wd_sz_list(iter_wd), wd_step_list(iter_wd), gf_list(iter_wd));
                else
                    tmp_corr_str = obj.compute_moving_trace_correlation(...
                        trace_1, trace_2, wd_sz_list(iter_wd), wd_step_list(iter_wd));
                end
                
                rho_mat(tmp_corr_str.wd_mid, iter_wd) = tmp_corr_str.rho;
                rho_p_mat(tmp_corr_str.wd_mid, iter_wd) = tmp_corr_str.rho_p;
                
                rho_d_mat(tmp_corr_str.wd_mid, iter_wd) = tmp_corr_str.rho_d;
                rho_d_p_mat(tmp_corr_str.wd_mid, iter_wd) = tmp_corr_str.rho_p_d;
                
                k_mat(tmp_corr_str.wd_mid, iter_wd) = tmp_corr_str.linear_regression_coeff(:, 1);
                b_mat(tmp_corr_str.wd_mid, iter_wd) = tmp_corr_str.linear_regression_coeff(:, 2);
            end
            mt_corr_str.raw.rho = rho_mat;
            mt_corr_str.raw.rho_p = rho_p_mat;
            mt_corr_str.raw.rho_d = rho_d_mat;
            mt_corr_str.raw.rho_d_p = rho_d_p_mat;
            mt_corr_str.raw.k = k_mat;
            mt_corr_str.raw.b = b_mat;
            fn = fieldnames(mt_corr_str.raw);
            for iter_fn = 1 : numel(fn)
                tmp_fn = fn{iter_fn};
                tmp_data = mt_corr_str.raw.(tmp_fn);
                mt_corr_str.(tmp_fn) = fillmissing(tmp_data, 'linear', 1, ...
                    'EndValues', 'none');
            end
        end
        
        function corr_str = compute_moving_trace_correlation(obj, trace_1, trace_2, wd_sz, wd_step, gf_sigma)
            if nargin < 6
                gf_sigma = [];
            end
            
            trace_len = numel(trace_1);
            assert(isvector(trace_1) && isvector(trace_2));
            assert(trace_len == numel(trace_2), 'Two input traces have different length');
            
            wd_start_ind = 1 : wd_step : (trace_len - wd_sz + 1);
            wd_end_ind = min(wd_start_ind + wd_sz - 1, trace_len);
            num_wds = numel(wd_start_ind);
            [rho, rho_p, rho_d, rho_d_p] = deal(nan(num_wds, 1));
            lr_coeff = nan(2, num_wds);
            % Apply gaussian filter?
            if ~isempty(gf_sigma)
                trace_1 = imgaussfilt(trace_1, gf_sigma);
                trace_2 = imgaussfilt(trace_2, gf_sigma);
            end
            
            for iter_wd = 1 : num_wds
                tmp_trace_1 = trace_1(wd_start_ind(iter_wd) : wd_end_ind(iter_wd));
                tmp_trace_2 = trace_2(wd_start_ind(iter_wd) : wd_end_ind(iter_wd));
                [rho(iter_wd), rho_p(iter_wd), rho_d(iter_wd), rho_d_p(iter_wd)] ...
                    = obj.compute_trace_correlation(...
                    tmp_trace_1, tmp_trace_2);
                lr_coeff(:, iter_wd) = obj.compute_linear_regression_parameter(tmp_trace_1, tmp_trace_2);
            end
            
            corr_str.wd_sz = wd_sz;
            corr_str.wd_step = wd_step;
            corr_str.wd_start = wd_start_ind.';
            corr_str.wd_end = wd_end_ind.';
            corr_str.wd_mid = round((wd_start_ind + wd_end_ind) / 2).';
            corr_str.rho = rho;
            corr_str.rho_p = rho_p;
            corr_str.rho_d = rho_d;
            corr_str.rho_p_d = rho_d_p;
            corr_str.linear_regression_coeff = lr_coeff.';
        end
        
    end
    
    methods(Static)
        function [rho_trace, p, rho_dtrace, d_p] = compute_trace_correlation(trace_1, trace_2)
            is_valid_Q = isfinite(trace_1) & isfinite(trace_2);
            if ~all(is_valid_Q)
                fprintf('%.2f of the trace is invalid\n', 1 - mean(is_valid_Q));
                trace_1 = trace_1(is_valid_Q);
                trace_2 = trace_2(is_valid_Q);
            end
            num_points = numel(trace_1);
            if num_points >= 3
                d_trace_1 = diff(trace_1);
                d_trace_2 = diff(trace_2);
                
                
                [rho_trace, p] = corrcoef(trace_1, trace_2);
                rho_trace = rho_trace(1,2);
                p = p(1,2);
                
                [rho_dtrace, d_p] = corrcoef(d_trace_1, d_trace_2);
                rho_dtrace = rho_dtrace(1,2);
                d_p = d_p(1,2);
            else
                rho_trace = nan;
                rho_dtrace = nan;
                p = nan;
                d_p = nan;
            end
        end
        
        function b = compute_linear_regression_parameter(x, y)
            if isrow(y)
                y = y.';
            end
            [num_data, ~] = size(x);
            assert(num_data == numel(y), 'x and y should have the same number of rows');
            x = cat(2, x, ones(num_data, 1));
            b = (x' * x) \ (x') * y;
        end
        
        function fig_hdl = vis_multiscale_moving_correlation(trace_1, trace_2, mt_corr_str, dt_s, ...
                trace_name)
            num_points = numel(trace_1);
            assert(num_points == numel(trace_2));
            plot_t = (1 : num_points) * dt_s;
            wd_list_s = mt_corr_str.window_s_list;
            
            % Visualization
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 2.5];
            num_subplot = 4;
            ax_hdl_0 = subplot(num_subplot, 1, 1);
            plot(ax_hdl_0, plot_t, trace_1);
            hold(ax_hdl_0, 'on');
            plot(ax_hdl_0, plot_t, trace_2);
            grid(ax_hdl_0, 'on');
            ax_hdl_0.YLabel.String = '\DeltaF/F';
            plot_t_idx = round(ax_hdl_0.XTick / dt_s);
            ax_hdl_0.XLabel.String = 'Time (s)';
            ax_hdl_0.YLim = prctile(cat(1, trace_1, trace_2), [0.1, 99.9]);
            
            if ~isempty(trace_name)
                legend(ax_hdl_0, trace_name);
            end
            
            ax_hdl_1 = subplot(num_subplot,1,2);
            im_1 = imagesc(ax_hdl_1, mt_corr_str.rho.');
            caxis(ax_hdl_1, [-1,1]);
            %         cbar_hdl_1 = colorbar(ax_hdl_1);
            %         cbar_hdl_1.Location = 'northoutside';
            %         cbar_hdl_1.Limits = [-1, 1];
            colormap(ax_hdl_1, 'jet');
            ax_hdl_1.Title.String = 'Trace correlation';
            
            ax_hdl_2 = subplot(num_subplot,1,3);
            im_2 = imagesc(ax_hdl_2, mt_corr_str.rho_d.');
            colormap(ax_hdl_2, 'jet');
            cbar_hdl_2 = colorbar(ax_hdl_2);
            cbar_hdl_2.Location = 'southoutside';
            caxis(ax_hdl_2, [-1, 1]);
            ax_hdl_2.Title.String = 'Derivative correlation';
            
            ax_hdl_3 = subplot(num_subplot, 1, 4);
            im_3 = imagesc(ax_hdl_3, mt_corr_str.k.');
            colormap(ax_hdl_3, 'jet');
            ax_hdl_3.Title.String = 'Slope';
            %             caxis(ax_hdl_3, prctile(mt_corr_str.k(:), [0.1, 99.9]));
            switch trace_name{2}
                case {'flow', 'Flow'}
                    caxis(ax_hdl_3, [-2, 2]);
                case {'flux', 'Flux'}
                    caxis(ax_hdl_3, [-8, 8]);
                otherwise
                    error('Unrecognized');
            end
            
            cbar_hdl_3 = colorbar(ax_hdl_3, 'Location', 'southoutside');
            
            [ax_hdl_1.YTick, ax_hdl_2.YTick, ax_hdl_3.YTick] = deal(1 : mt_corr_str.num_windows);
            [ax_hdl_1.YTickLabel, ax_hdl_2.YTickLabel, ax_hdl_3.YTickLabel] = ...
                deal(arrayfun(@(x) num2str(x, '%d'), wd_list_s, 'UniformOutput', false));
            [ax_hdl_1.YLabel.String, ax_hdl_2.YLabel.String, ax_hdl_3.YLabel.String] = ...
                deal('Window size (s)');
            [ax_hdl_1.XTick, ax_hdl_2.XTick, ax_hdl_3.XTick] = deal(plot_t_idx);
            [ax_hdl_1.XTickLabel, ax_hdl_2.XTickLabel, ax_hdl_3.XTickLabel] = deal([]);
            
            %             k_limit = max(abs(tmp_wds_str.k(:)));
            %             cbar_hdl_3.Limits = [-k_limit, k_limit];
            
        end
    end
    %% Analyze fluctuation near the heart beat frequency 
    methods
        
    end 
    
    methods(Static)
        function stat = compute_trace_modulation_in_frequency_range(data, f_range_Hz, f_sample_Hz, num_drop_tps)
            if nargin < 4
                num_drop_tps = 0;
            end
            f_range_Hz = min(f_range_Hz, f_sample_Hz);
            r_or_bp = bandpass(data, f_range_Hz, f_sample_Hz);
            r_or_bp = r_or_bp((1 + num_drop_tps) : (end - num_drop_tps));
            stat = fun_analysis_get_basic_statistics(abs(r_or_bp));
        end        
    end
    %% Visualization
    methods        
        % Coherence
        function fig_hdl = vis_pairwise_coherence(obj, save_fig_Q)
            if nargin < 2
                save_fig_Q = true;
            end
            name_list = {'Radius', 'Flux', 'Speed'};
            num_pair = numel(name_list);
            vis_f_list = obj.frequency.r.psd.f;
            
            for vis_idx_1 = 1 : (num_pair - 1)
                for vis_idx_2 = (vis_idx_1 + 1) : num_pair
                    
                    vis_name_1 = name_list{vis_idx_1};
                    vis_name_2 = name_list{vis_idx_2};
                    
                    vis_abs = obj.frequency.coherence.abs{vis_idx_1,vis_idx_2};
                    vis_ang = obj.frequency.coherence.arg{vis_idx_1,vis_idx_2};
                    
                    fig_hdl = figure;
                    ax_hdl = axes(fig_hdl);
                    yyaxis(ax_hdl, 'left');
                    plot(ax_hdl, vis_f_list, vis_abs, 'LineWidth', 1);
                    ax_hdl.XScale = deal('log');
                    %     ax_hdl.YScale = deal('log');
                    ax_hdl.XLabel.String = deal('Frequency (Hz)');
                    ax_hdl.YLabel.String = deal('|Coherence|');
                    ax_hdl.XLim = [obj.frequency.r.psd.f_resolve_Hz, vis_f_list(end)];
                    yyaxis(ax_hdl, 'right');
                    plot(ax_hdl, vis_f_list, vis_ang * 180 / pi, 'LineWidth', 1);
                    ax_hdl.YLabel.String = 'Phase (degrees)';
                    grid(ax_hdl, 'on');
                    ax_hdl.Title.String = sprintf('%s vs %s', vis_name_1, vis_name_2);
                    
                    if save_fig_Q
                        fig_fp = fullfile(obj.vis_root_folder, 'average', 'Coherence', ...
                            sprintf('%s_coherence_%s_vs_%s.png', obj.vis_fp_template, vis_name_1, vis_name_2));
                        fun_print_image_in_several_formats(fig_hdl, fig_fp);
                        delete(fig_hdl);
                    end
                end
            end
        end
    end
    
    
    methods(Static)
        
        % Power spectrum density
        function fig_hdl = vis_power_spectrum_densities(psd_cell, psd_label_cell)
            vis_cell = cat(1, psd_cell, psd_label_cell);
            num_vis = size(vis_cell, 2);
            
            fig_hdl = figure;
            fig_hdl.Position(4) = fig_hdl.Position(4) * 2;
            for iter_vis = 1 : num_vis
                ax_hdl = subplot(num_vis, 1, iter_vis);
                tmp_psd = vis_cell{1, iter_vis};
                tmp_title = vis_cell{2, iter_vis};
                
                plt_hdl = plot(ax_hdl, tmp_psd.f, tmp_psd.psd, 'LineWidth', 1);
                ax_hdl.XScale = deal('log');
                ax_hdl.YScale = deal('log');
                if iter_vis == num_vis
                    ax_hdl.XLabel.String = deal('Frequency (Hz)');
                end
                ax_hdl.YLabel.String = deal('Spectral Power Density');
                ax_hdl.XLim = [tmp_psd.f_resolve_Hz, tmp_psd.f(end)];
                grid(ax_hdl, 'on');
                tmp_cpsd = cumsum(tmp_psd.psd);
                tmp_cpsd = tmp_cpsd ./ tmp_cpsd(end);
                hold(ax_hdl, 'on');
                yyaxis(ax_hdl, 'right');
                plot(ax_hdl, tmp_psd.f, tmp_cpsd);
                ax_hdl.YLabel.String = 'Normalized cumulative power';
                ax_hdl.YLim(1) = 0;
                ax_hdl.Title.String = tmp_title;
            end
        end
        
        function  fig_hdl = vis_average_traces(r_trace, flux, t_trace, ...
                flow_mov_wd_sz)
            if nargin < 4
                flow_mov_wd_sz = [];
            end
            
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 2];
            ax_hdl_1 = subplot(3,1,1);
            ax_hdl_2 = subplot(3,1,2);
            ax_hdl_3 = subplot(3,1,3);
            
            plot(ax_hdl_1, t_trace, r_trace);
            ax_hdl_1.YLabel.String = 'Radius (\mum)';
            ax_hdl_1.YLim(1) = 0;
            grid(ax_hdl_1, 'on');
            hold(ax_hdl_1, 'on');
            if ~isempty(flow_mov_wd_sz)
                r_sm = movmean(r_trace, flow_mov_wd_sz, 'omitnan');
                plot(ax_hdl_1, t_trace, r_sm);
            end
            
            plot(ax_hdl_2, t_trace, flux);
            ax_hdl_2.YLabel.String = 'RBC Flux (RBC/s)';
            ax_hdl_2.YLim(1) = 0;
            grid(ax_hdl_2, 'on');
            if ~isempty(flow_mov_wd_sz)
                hold(ax_hdl_2, 'on');
                flux_sm = movmean(flux, flow_mov_wd_sz, 'omitnan');
                plot(ax_hdl_2, t_trace, flux_sm);
            end
            
            yyaxis(ax_hdl_3, 'left');
            if ~isempty(flow_mov_wd_sz)
                plot(ax_hdl_3, t_trace, r_sm);
            else
                plot(ax_hdl_3, t_trace, r_trace);
            end
            ax_hdl_3.YLabel.String = 'Radius (\mum)';
            %                 ax_hdl_4.YLim(1) = 0;
            yyaxis(ax_hdl_3, 'right');
            if ~isempty(flow_mov_wd_sz)
                plot(ax_hdl_3, t_trace, flux_sm);
            else
                plot(ax_hdl_3, t_trace, flux);
            end
            ax_hdl_3.YLabel.String = 'Flux (RBC/s)';
            %                 ax_hdl_4.YLim(1) = 0;
            ax_hdl_3.XLabel.String = 'Time (s)';
            grid(ax_hdl_3, 'on');
            %                 if ~isempty(flow_mov_wd_sz)
            %                     hold(ax_hdl_4, 'on');
            %                     plot(ax_hdl_4, flux_sm, 'LineWidth', 1);
            %                 end
        end
        
        function fig_hdl = vis_radius_and_flux_with_r_phase(t_trace, r_trace, ...
                flux_trace, phase_trace)
            fig_hdl = figure;
            ax_hdl_1 = subplot(2,1,1);
            yyaxis(ax_hdl_1, 'left');
            plot(ax_hdl_1, t_trace, r_trace);
            ax_hdl_1.YLabel.String = 'Radius (\mum)';
            hold(ax_hdl_1, 'on');
            yyaxis(ax_hdl_1, 'right');
            plot(ax_hdl_1, t_trace, phase_trace);
            grid(ax_hdl_1, 'on');
            ax_hdl_1.YLabel.String = 'Phase (rad)';
            ax_hdl_1.XLabel.String = 'Time (s)';
            
            ax_hdl_2 = subplot(2,1,2);
            yyaxis(ax_hdl_2, 'left');
            plot(ax_hdl_2, t_trace, flux_trace);
            ax_hdl_2.YLabel.String = 'Flux (RBC/s)';
            hold(ax_hdl_2, 'on');
            yyaxis(ax_hdl_2, 'right');
            plot(ax_hdl_2, t_trace, phase_trace);
            grid(ax_hdl_2, 'on');
            ax_hdl_2.YLabel.String = 'Phase (rad)';
            ax_hdl_2.XLabel.String = 'Time (s)';
        end
        
        function fig_hdl = vis_radius_vs_speed(t_trace, r_trace, y_trace, ...
                x_name_1, y_name_1, x_name_2, y_name_2)
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,1];
            
            ax_hdl = subplot(1,2,1);
            yyaxis(ax_hdl, 'left');
            plot(ax_hdl, t_trace, r_trace);
            ax_hdl.YLabel.String = x_name_1;
            hold(ax_hdl, 'on');
            ax_hdl.XLabel.String = 'Time (s)';
            yyaxis(ax_hdl, 'right');
            plot(ax_hdl, t_trace, y_trace);
            ax_hdl.YLabel.String = y_name_1;
            ax_hdl_2 = subplot(1,2,2);
            grid(ax_hdl, 'on');
            
            r_mean = mean(r_trace);
            v_mean = mean(y_trace);
            
            r_sm_n = (r_trace - r_mean) / r_mean;
            v_sm_n = (y_trace - v_mean) / v_mean;
            dvdr = fitlm(r_sm_n, v_sm_n, 'Intercept', false);
            histogram2(ax_hdl_2, r_sm_n, v_sm_n, 24, ...
                'DisplayStyle', 'tile');
            ax_hdl_2.XLabel.Interpreter = 'latex';
            ax_hdl_2.YLabel.Interpreter = 'latex';
            ax_hdl_2.XLabel.String = x_name_2;
            ax_hdl_2.YLabel.String = y_name_2;
            ax_hdl_2.XLabel.FontSize = 14;
            ax_hdl_2.YLabel.FontSize = 14;
            ax_hdl_2.XLim = [-1, 1] .* max(abs(ax_hdl_2.XLim));
            ax_hdl_2.YLim = [-1, 1] .* max(abs(ax_hdl_2.YLim));
            hold(ax_hdl_2, 'on');
            plt_x = linspace(ax_hdl_2.XLim(1), ax_hdl_2.XLim(2), 100);
            plt_y = dvdr.Coefficients.Estimate * plt_x;
            line_hdl = plot(ax_hdl_2, plt_x, plt_y, 'LineWidth', 2);
            legend(ax_hdl_2, line_hdl, sprintf('Slope: %.3f\nR^2: %.3f',...
                dvdr.Coefficients.Estimate, dvdr.Rsquared.Adjusted));
        end
    end
end