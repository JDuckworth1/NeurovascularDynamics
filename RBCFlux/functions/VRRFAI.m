classdef VRRFAI < handle
    % ROI integration
    properties
        raw_data VRRFA
        data struct
        dataset_name char
        stack char
        roi_idx double
        record_idx double
        vessel_type char
        title_prefix char
        parameter struct
        
        radius struct
        flux struct
        speed struct        
        frequency struct
        
        flux_vs_phase struct
        radius_vs_phase struct
        radius_vs_flux struct
        radius_vs_speed struct
        radius_vs_flux_intercept struct
        
        flux_vs_flow struct
        
        vis_root_folder char
        filepath char
        vis_fp_template char
        
        norm_r double
        norm_flux double
        r_trace double
        flux_trace double
    end
    
    methods
        function obj = VRRFAI(vrrfa_array, vsl_type)
            obj.raw_data = vrrfa_array;
            obj.parameter = cat(1, obj.raw_data.parameters);
            
            obj.stack = unique(cat(1, vrrfa_array.stack), 'rows');
            obj.dataset_name = unique(cat(1, vrrfa_array.dataset_name), 'rows');
            obj.roi_idx = unique([vrrfa_array.roi_idx]);
            obj.record_idx = [vrrfa_array.record_idx];
            obj.vessel_type = vsl_type;
            
            DataManager = FileManager;
            obj.vis_root_folder = fullfile(sprintf('%s_sl', DataManager.fp_visualization_folder(...
                obj.dataset_name, obj.stack)), sprintf('ROI_%03d', obj.roi_idx));
            obj.vis_fp_template = sprintf('%s_%s_ROI_%03d', obj.dataset_name, ...
                obj.stack, obj.roi_idx);
            obj.title_prefix = sprintf('ROI %03d %s Vessel', obj.roi_idx, ...
                obj.vessel_type);
        end
        
        function parse_parameters(obj, ppp_str)
            obj.parameter = ppp_str;
            obj.parameter.num_phase_edge = obj.parameter.num_phase_bin + 1;
            obj.parameter.num_r_edge = obj.parameter.num_r_bin + 1;
            obj.parameter.num_flux_edge = obj.parameter.num_flux_bin + 1;
        end
        
        function obj = merge_data(obj)
            fn_list = {'flux.per_frame.data', 'flux.per_frame_mm.data', ...
                'flux.per_frame_mm.vm_int', 'flux.per_frame_mm.vm_int_stat', ...
                'speed.per_frame.data', 'speed.per_frame_mm.data', ...
                'speed.per_frame_mm.vm_int', 'speed.per_frame_mm.vm_int_stat', ...
                'speed.per_frame_mm.est_flux', 'speed.per_frame_mm.est_flux_n', ...
                'radius.r_or.data', 'radius.r_or.phase', ...
                'radius.r_or_mm.data', 'radius.r_or_mm.phase', 'radius.r_or_mm.phase_l'};
            new_str = struct;
            for i_fn = 1 : numel(fn_list)
                tmp_fn = fn_list{i_fn};
                tmp_data = arrayfun(@(x) fun_getfield(x, tmp_fn), obj.raw_data,...
                    'UniformOutput', false);
                tmp_data = cat(1, tmp_data{:});
                tmp_new_fn = strrep(tmp_fn, '.', '_');
                new_str.(tmp_new_fn) = tmp_data;
            end
            obj.data = new_str;
        end
        
        function obj = merge_radius_data(obj)
            % Get data
            obj.radius = struct;
            fn_list = {'r_or', 'r_or_mm'};
            for i = 1 : numel(fn_list)
                tmp_fn = fn_list{i};
                tmp_str_cell = arrayfun(@(x) x.radius.(tmp_fn), obj.raw_data, 'UniformOutput', false);
                obj.radius.(tmp_fn) = obj.merge_data_str(tmp_str_cell);
            end            
        end
        
        function obj = merge_flux_and_speed_data(obj)
            fn_list_1 = {'speed', 'flux'};
            fn_list_2 = {'per_frame', 'per_frame_mm'};
            for i = 1 : numel(fn_list_1)
                tmp_fn1 = fn_list_1{i};
                obj.(tmp_fn1) = struct;
                for j = 1 : numel(fn_list_2)
                    tmp_fn2 = fn_list_2{j};
                    tmp_str_cell = arrayfun(@(x) x.(tmp_fn1).(tmp_fn2), obj.raw_data, ...
                        'UniformOutput', false);
                    obj.(tmp_fn1).(tmp_fn2) = obj.merge_data_str(...
                        tmp_str_cell);
                end
            end
        end
        
        
        function flux_in_phase_bin_str = analyze_vs_phase(obj, phase_trace, flux_trace)
            flux_in_phase_bin_str = fun_analysis_get_y_stat_in_x_bin(...
                phase_trace, flux_trace, obj.parameter.phase_bin_edge);
            tmp_fit_y = flux_in_phase_bin_str.y_median;
            tmp_fit_x = cos(flux_in_phase_bin_str.x_bin_val);
            tmp_fit = fitlm(tmp_fit_x, tmp_fit_y);
            flux_in_phase_bin_str.fit_med.c = tmp_fit.Coefficients.Estimate(1);
            flux_in_phase_bin_str.fit_med.A = tmp_fit.Coefficients.Estimate(2);
            flux_in_phase_bin_str.fit_med.R2 = tmp_fit.Rsquared.Adjusted;
            % Modulation depth
            flux_in_phase_bin_str.fit_med.A2c = flux_in_phase_bin_str.fit_med.A / ...
                flux_in_phase_bin_str.fit_med.c;
            
            is_valid_Q = isfinite(phase_trace) & isfinite(flux_trace);
            tmp_fit_raw = fitlm(cos(phase_trace(is_valid_Q)), flux_trace(is_valid_Q));
            flux_in_phase_bin_str.fit_raw.c = tmp_fit_raw.Coefficients.Estimate(1);
            flux_in_phase_bin_str.fit_raw.A = tmp_fit_raw.Coefficients.Estimate(2);
            flux_in_phase_bin_str.fit_raw.R2 = tmp_fit_raw.Rsquared.Adjusted;
            % Modulation depth
            flux_in_phase_bin_str.fit_raw.A2c = flux_in_phase_bin_str.fit_raw.A / ...
                flux_in_phase_bin_str.fit_raw.c;
            
        end
        
        function y_in_x_bin_str = analyze_x_vs_y(obj, x_trace, y_trace)
            x_stat = fun_analysis_get_basic_statistics(x_trace);
            x_edge = linspace(x_stat.prctile_val(2), x_stat.prctile_val(end-1),...
                obj.parameter.num_r_edge);
            y_in_x_bin_str = fun_analysis_get_y_stat_in_x_bin(...
                x_trace, y_trace, x_edge);
            tmp_fit_y = y_in_x_bin_str.y_median;
            tmp_fit_x = y_in_x_bin_str.x_bin_val;
            tmp_fit = fitlm(tmp_fit_x, tmp_fit_y);
            y_in_x_bin_str.fit_med.b = tmp_fit.Coefficients.Estimate(1);
            y_in_x_bin_str.fit_med.k = tmp_fit.Coefficients.Estimate(2);
            y_in_x_bin_str.fit_med.R2 = tmp_fit.Rsquared.Adjusted;
            
            tmp_fit_all = fitlm(x_trace, y_trace);
            y_in_x_bin_str.fit_raw.b = tmp_fit_all.Coefficients.Estimate(1);
            y_in_x_bin_str.fit_raw.k = tmp_fit_all.Coefficients.Estimate(2);
            y_in_x_bin_str.fit_raw.R2 = tmp_fit_all.Rsquared.Adjusted;
        end
        
        function y_in_x_bin_str = analyze_x_vs_y_smrange(obj, x_trace, y_trace)
            x_stat = fun_analysis_get_basic_statistics(x_trace);
            x_edge = linspace(x_stat.prctile_val(2), x_stat.prctile_val(end-1),...
                obj.parameter.num_r_edge);
            y_in_x_bin_str = fun_analysis_get_y_stat_in_x_bin(...
                x_trace, y_trace, x_edge);
            tmp_fit_y = y_in_x_bin_str.y_median;
            tmp_fit_x = y_in_x_bin_str.x_bin_val;
            tmp_fit = fitlm(tmp_fit_x, tmp_fit_y);
            y_in_x_bin_str.fit_med.b = tmp_fit.Coefficients.Estimate(1);
            y_in_x_bin_str.fit_med.k = tmp_fit.Coefficients.Estimate(2);
            y_in_x_bin_str.fit_med.R2 = tmp_fit.Rsquared.Adjusted;
            
            tmp_fit_all = fitlm(x_trace, y_trace);
            y_in_x_bin_str.fit_raw.b = tmp_fit_all.Coefficients.Estimate(1);
            y_in_x_bin_str.fit_raw.k = tmp_fit_all.Coefficients.Estimate(2);
            y_in_x_bin_str.fit_raw.R2 = tmp_fit_all.Rsquared.Adjusted;
        end
        
        function dyn_in_dxn_bin_str = analyze_normalized_changes(obj, dx_n_trace, dy_n_trace)
            x_stat = fun_analysis_get_basic_statistics(dx_n_trace);
            x_edge = linspace(x_stat.prctile_val(2), x_stat.prctile_val(end-1),...
                obj.parameter.num_r_edge);
            dyn_in_dxn_bin_str = fun_analysis_get_y_stat_in_x_bin(...
                dx_n_trace, dy_n_trace, x_edge);
            tmp_fit_y = dyn_in_dxn_bin_str.y_median; %Median of y in every x bin
            tmp_fit_x = dyn_in_dxn_bin_str.x_bin_val; %X bin values
            tmp_fit = fitlm(tmp_fit_x, tmp_fit_y, 'Intercept', false);
            dyn_in_dxn_bin_str.fit_med.k = tmp_fit.Coefficients.Estimate(1);
            dyn_in_dxn_bin_str.fit_med.R2 = tmp_fit.Rsquared.Adjusted;
            
            tmp_fit_all = fitlm(dx_n_trace, dy_n_trace, 'Intercept', false);
            dyn_in_dxn_bin_str.fit_raw.k = tmp_fit_all.Coefficients.Estimate(1);
            dyn_in_dxn_bin_str.fit_raw.R2 = tmp_fit_all.Rsquared.Adjusted;
        end
        
        function dyn_in_dxn_bin_str = analyze_normalized_changes_intercept(obj, dx_n_trace, dy_n_trace)
            x_stat = fun_analysis_get_basic_statistics(dx_n_trace);
            x_edge = linspace(x_stat.prctile_val(2), x_stat.prctile_val(end-1),...
                obj.parameter.num_r_edge);
            dyn_in_dxn_bin_str = fun_analysis_get_y_stat_in_x_bin(...
                dx_n_trace, dy_n_trace, x_edge);
            tmp_fit_y = dyn_in_dxn_bin_str.y_median; %Median of y in every x bin
            tmp_fit_x = dyn_in_dxn_bin_str.x_bin_val; %X bin values
            tmp_fit = fitlm(tmp_fit_x, tmp_fit_y, 'Intercept', true);
            dyn_in_dxn_bin_str.fit_med.int = tmp_fit.Coefficients.Estimate(1);
            dyn_in_dxn_bin_str.fit_med.k = tmp_fit.Coefficients.Estimate(2);
            dyn_in_dxn_bin_str.fit_med.R2 = tmp_fit.Rsquared.Adjusted;
            dyn_in_dxn_bin_str.fit_med.k_SE = tmp_fit.Coefficients.SE(2);
            dyn_in_dxn_bin_str.fit_med.int_SE = tmp_fit.Coefficients.SE(1);           
            
            tmp_fit_all = fitlm(dx_n_trace, dy_n_trace, 'Intercept', true);
            dyn_in_dxn_bin_str.fit_raw.int = tmp_fit_all.Coefficients.Estimate(1);
            dyn_in_dxn_bin_str.fit_raw.k = tmp_fit_all.Coefficients.Estimate(2);
            dyn_in_dxn_bin_str.fit_raw.R2 = tmp_fit_all.Rsquared.Adjusted;
            
            dyn_in_dxn_bin_str.fit_raw.k_SE = tmp_fit_all.Coefficients.SE(2);
            dyn_in_dxn_bin_str.fit_raw.int_SE = tmp_fit_all.Coefficients.SE(1);
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
        
        
                
    end
    %% Frequency
    methods
        function result = compute_fractional_power_in_frequency_range(obj, f_range_Hz)
            % Overestimate
            num_raw_data = numel(obj.raw_data);
            analyze_fields = {'r', 'flux', 'speed'};
            result = struct;
            for i = 1 : num_raw_data
                tmp_hdl = obj.raw_data(i).frequency;
                for j = 1 : numel(analyze_fields)
                    tmp_field = analyze_fields{j};
                    if i == 1
                        result.(tmp_field).data = zeros(num_raw_data, 1);                        
                    end                    
                    tmp_psd = tmp_hdl.(tmp_field).psd;
                    tmp_f_in_range_Q = tmp_psd.f >= f_range_Hz(1) & ...
                        tmp_psd.f <= f_range_Hz(2);
                    result.(tmp_field).data(i) = sum(tmp_psd.psd(tmp_f_in_range_Q)) ./ ...
                        sum(tmp_psd.psd);                    
                    if i == num_raw_data
                        result.(tmp_field).stat = fun_analysis_get_basic_statistics(result.(tmp_field).data);                        
                    end
                end                
            end                        
        end
        
        function obj = compute_radius_modulation_in_frequency_range(obj, f_range_Hz, num_drop_tps)
            if nargin < 3
                num_drop_tps = 0;
            end
            num_raw_data = numel(obj.raw_data);
            stat = struct;
            [stat.data_mean, stat.data_median] = deal(zeros(num_raw_data, 1));
            for i = 1 : num_raw_data
                tmp_hdl = obj.raw_data(i);
                tmp_hdl.radius.r_or.stat_in_fband = tmp_hdl.compute_trace_modulation_in_frequency_range(...
                    tmp_hdl.radius.r_or.data, f_range_Hz, tmp_hdl.parameters.frame_rate_Hz, num_drop_tps);

                stat.data_mean(i) = tmp_hdl.radius.r_or.stat_in_fband.mean / ...
                    tmp_hdl.radius.r_or.stat.mean;
                stat.data_median(i) =  tmp_hdl.radius.r_or.stat_in_fband.median / ...
                    tmp_hdl.radius.r_or.stat.median;
            end
            stat.stat_mean = fun_analysis_get_basic_statistics(stat.data_mean);
            stat.stat_median = fun_analysis_get_basic_statistics(stat.data_median);      
            obj.radius.r_or.rm_in_fband = stat;
        end
    end    
    %%    
    methods(Static)
        function result_str = merge_data_str(radius_str_cell)
            radius_str_cell = cat(1, radius_str_cell{:});
            fn = fieldnames(radius_str_cell);
            for i = 1 : numel(fn)
                tmp_fn = fn{i};
                if isnumeric(radius_str_cell(1).(tmp_fn)) || strcmpi(tmp_fn, 'vm_int') ...
                        || istable(radius_str_cell(1).(tmp_fn))
                   tmp_data = cat(1, radius_str_cell.(tmp_fn)); 
                end
                result_str.(tmp_fn) = tmp_data;
            end
            result_str.stat = fun_analysis_get_basic_statistics(...
                result_str.data);
            result_str.data_n = (result_str.data - result_str.stat.mean) ./ ...
                result_str.stat.mean;
            result_str.stat_n = fun_analysis_get_basic_statistics(result_str.data_n);
        end
        
        function vm_int_stat = select_vasomition_interval_stat(vm_int_stat, ...
                min_duration_frac)
            is_valid_Q = vm_int_stat.rf_frac(:, 1) > min_duration_frac &...
                vm_int_stat.rf_frac(:, 1) < 1 - min_duration_frac;
            vm_int_stat = vm_int_stat(is_valid_Q, :);
        end
    end    
    %% Visualization
    methods(Static)
        function fig_hdl = vis_vasomotion_interval_stat(vm_int_stat)           
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 2;
            ax_hdl_1 = subplot(2,3,1);
            scatter(ax_hdl_1, vm_int_stat.dfdr_n_k, vm_int_stat.r2f_r2);
            grid(ax_hdl_1, 'on');
            ax_hdl_1.XLabel.String = 'Slope';
            ax_hdl_1.YLabel.String = 'R^2';
            ax_hdl_2 = subplot(2,3,2);
            scatter(ax_hdl_2, vm_int_stat.r_range, vm_int_stat.dfdr_n_k);
            grid(ax_hdl_2, 'on');
            ax_hdl_2.XLabel.String = '\Delta r (\mum)';
            ax_hdl_2.YLabel.String = 'Slope';
            ax_hdl_3 = subplot(2,3,3);
            scatter(ax_hdl_3, vm_int_stat.r_range, vm_int_stat.r2f_r2);
            grid(ax_hdl_3, 'on');
            ax_hdl_3.XLabel.String = '\Delta r (\mum)';
            ax_hdl_3.YLabel.String = 'R^2';
            ax_hdl_4 = subplot(2,3,4);
            scatter(ax_hdl_4, vm_int_stat.r_range, vm_int_stat.durations);
            grid(ax_hdl_4, 'on');
            ax_hdl_4.XLabel.String = '\Deltar(\mum)';
            ax_hdl_4.YLabel.String = 'Duration (s)';
            ax_hdl_5 = subplot(2,3,5);
            scatter(ax_hdl_5, vm_int_stat.r_range ./ vm_int_stat.r_mean, vm_int_stat.r2f_r2);
            grid(ax_hdl_5, 'on');
            ax_hdl_5.XLabel.String = '\Deltar/\langler\rangle';
            ax_hdl_5.YLabel.String = 'R^2';            
            ax_hdl_6 = subplot(2,3,6);
            scatter(ax_hdl_6, vm_int_stat.r_range ./ vm_int_stat.r_mean, vm_int_stat.dfdr_n_k);
            grid(ax_hdl_6, 'on');
            ax_hdl_6.XLabel.String = '\Deltar/\langler\rangle';
            ax_hdl_6.YLabel.String = 'Slope';
        end
        
        function fig_hdl = vis_vasomotion_interval_stat_hist2(vm_int_stat)
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 2;
            ax_hdl_1 = subplot(2,3,1);
            histogram2(ax_hdl_1, vm_int_stat.dfdr_n_k, vm_int_stat.r2f_r2, ...
                'DisplayStyle', 'tile');
            grid(ax_hdl_1, 'on');
            ax_hdl_1.XLabel.String = 'Slope';
            ax_hdl_1.YLabel.String = 'R^2';
            ax_hdl_2 = subplot(2,3,2);
            histogram2(ax_hdl_2, vm_int_stat.r_range, vm_int_stat.dfdr_n_k, ...
                'DisplayStyle', 'tile');
            grid(ax_hdl_2, 'on');
            ax_hdl_2.XLabel.String = '\Delta r (\mum)';
            ax_hdl_2.YLabel.String = 'Slope';
            ax_hdl_3 = subplot(2,3,3);
            histogram2(ax_hdl_3, vm_int_stat.r_range, vm_int_stat.r2f_r2,...
                'DisplayStyle', 'tile');
            grid(ax_hdl_3, 'on');
            ax_hdl_3.XLabel.String = '\Delta r (\mum)';
            ax_hdl_3.YLabel.String = 'R^2';
            ax_hdl_4 = subplot(2,3,4);
            histogram2(ax_hdl_4, vm_int_stat.r_range, vm_int_stat.durations,...
                'DisplayStyle', 'tile');
            grid(ax_hdl_4, 'on');
            ax_hdl_4.XLabel.String = '\Deltar(\mum)';
            ax_hdl_4.YLabel.String = 'Duration (s)';
            ax_hdl_5 = subplot(2,3,5);
            histogram2(ax_hdl_5, vm_int_stat.r_range ./ vm_int_stat.r_mean, vm_int_stat.r2f_r2,...
                'DisplayStyle', 'tile');
            grid(ax_hdl_5, 'on');
            ax_hdl_5.XLabel.String = '\Deltar/\langler\rangle';
            ax_hdl_5.YLabel.String = 'R^2';            
            ax_hdl_6 = subplot(2,3,6);
            histogram2(ax_hdl_6, vm_int_stat.r_range ./ vm_int_stat.r_mean, vm_int_stat.dfdr_n_k,...
                'DisplayStyle', 'tile');
            grid(ax_hdl_6, 'on');
            ax_hdl_6.XLabel.String = '\Deltar/\langler\rangle';
            ax_hdl_6.YLabel.String = 'Slope';
        end
        
        function fig_hdl = vis_vasomotion_interval_stat_v2(vm_int_stat, plt_title)
            if nargin < 2
                plt_title = [];
            end
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            sc_hdl = scatter(ax_hdl, vm_int_stat.dfdr_n_k, vm_int_stat.r2f_r2, ...
                vm_int_stat.r_mean * 10, ...
                vm_int_stat.r_range ./ vm_int_stat.r_mean, 'filled');
            cbar_hdl = colorbar(ax_hdl);
            colormap(ax_hdl, 'jet');
            sc_hdl.AlphaData = 0.5 * ones(size(sc_hdl.XData));
            ax_hdl.XLabel.String = 'Slope';
            ax_hdl.YLabel.String = 'R^2';
            ax_hdl.YLim = [0,1];
            cbar_hdl.Label.String = '\Deltar/\langler\rangle';
            grid(ax_hdl, 'on');
            if ~isempty(plt_title)
                ax_hdl.Title.String = plt_title; 
            end
        end
        
        function [fig_hdl, varargout] = vis_xy_joint_dist_in_hist(x, y, x_edge, y_edge, ...
                binned_stat, x_name, y_name, plt_title)
            if nargin < 8
                plt_title = [];
            end
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            if ~isempty(x_edge) && ~isempty(y_edge)
                hist_hdl = histogram2(ax_hdl, x, y, x_edge, y_edge, 'DisplayStyle', 'tile');
            else
                hist_hdl = histogram2(ax_hdl, x, y, 'DisplayStyle', 'tile');
                x_edge = hist_hdl.XBinEdges;
%                 y_edge = hist_hdl.YBinEdges;
            end
            if isempty(binned_stat)
                binned_stat = fun_analysis_get_y_stat_in_x_bin(x, y, ...
                    x_edge);
            end
            hold(ax_hdl, 'on');
            eb_hdl = errorbar(ax_hdl, binned_stat.x_bin_val, binned_stat.y_median, ...
                binned_stat.y_median - binned_stat.y_prctile(:, 3), ...
                binned_stat.y_prctile(:, 5) - binned_stat.y_median, ...
                'LineWidth', 1.5);
            ax_hdl.XLabel.String = x_name;
            ax_hdl.YLabel.String = y_name;
            
            if ~isempty(plt_title)
                ax_hdl.Title.String = plt_title;
            end
            if nargout > 1
                varargout{1} = ax_hdl;
            end
        end
        
        function fig_hdl = vis_y_vs_phase_in_hist(x, y, x_edge, y_edge, ...
                binned_stat, x_name, y_name, plt_title)
            [fig_hdl, ax_hdl] = VRRFAI.vis_xy_joint_dist_in_hist(x, y, x_edge, ...
                y_edge, binned_stat, x_name, y_name, plt_title);
            hold(ax_hdl, 'on');
            fit_hdl = plot(ax_hdl, binned_stat.x_bin_val, ...
                binned_stat.fit_med.c + binned_stat.fit_med.A .* cos(binned_stat.x_bin_val), ...
                'k', 'LineWidth', 1);
            legend(ax_hdl, fit_hdl, sprintf('A cos(\\theta) + c\nA: %.2e\nc: %.2e\nA/c: %.2e\nR^2: %.2f', ...
                binned_stat.fit_med.A, binned_stat.fit_med.c, binned_stat.fit_med.A2c, ...
                binned_stat.fit_med.R2));
        end
        
        function fig_hdl = vis_radius_vs_flux(x, y, x_edge, y_edge, ...
                binned_stat, x_name, y_name, plt_title)
            [fig_hdl, ax_hdl] = VRRFAI.vis_xy_joint_dist_in_hist(x, y, x_edge, ...
                y_edge, binned_stat, x_name, y_name, plt_title);
            hold(ax_hdl, 'on');
            if isfield(binned_stat.fit_med, 'b')
                fit_hdl = plot(ax_hdl, binned_stat.x_bin_val, ...
                    binned_stat.fit_med.b + binned_stat.fit_med.k .* binned_stat.x_bin_val, ...
                    'k', 'LineWidth', 1);
                legend(ax_hdl, fit_hdl, sprintf('k r + b\nk %.2e\nb: %.2e\nR^2: %.2f', ...
                    binned_stat.fit_med.k, binned_stat.fit_med.b, binned_stat.fit_med.R2), ...
                    'Location', 'best');
            else
                fit_hdl = plot(ax_hdl, binned_stat.x_bin_val, ...
                     binned_stat.fit_med.k .* binned_stat.x_bin_val, ...
                    'k', 'LineWidth', 1);
                legend(ax_hdl, fit_hdl, sprintf('k r\nk %.2e\nR^2: %.2f', ...
                    binned_stat.fit_med.k, binned_stat.fit_med.R2), ...
                    'Location', 'best');
            end
        end
        
        function fig_hdl = vis_flux_vs_phase_controlled_by_radius()
            % Bin data by radius
            num_r_bin = 6;
            r_bin_edge = linspace(r_stat.prctile_val(2), r_stat.prctile_val(end-1), num_r_bin + 1);
            
            [idx_list, r_bin_val] = fun_bin_data_to_idx_list_by_edges(r_trace, r_bin_edge);
            
            binned_phase = fun_bin_data_to_cells_by_ind(phase_trace, idx_list);
            
            flux_lim = [flux_stat.min, flux_stat.max];
            binned_flux = fun_bin_data_to_cells_by_ind(flux_trace, idx_list);
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 3;
            %             ax_hdl = axes(fig_hdl);
            for i_bin = 1 : num_r_bin
                tmp_ax_hdl = subplot(2,3,i_bin);
                
                tmp_phase = binned_phase{i_bin};
                tmp_flux = binned_flux{i_bin};
                if ~isempty(tmp_phase)
                    tmp_is_valid_Q = isfinite(tmp_phase) & isfinite(tmp_flux);
                    tmp_phase = tmp_phase(tmp_is_valid_Q);
                    tmp_flux = tmp_flux(tmp_is_valid_Q);
                end
                
                if ~isempty(tmp_phase)
                    [tmp_phase_bin, tmp_phase_val] = fun_bin_data_to_idx_list_by_edges(tmp_phase, ...
                        phase_edge);
                    tmp_flux_bin_stat = fun_analysis_get_basic_statistics_in_bins(tmp_flux, ...
                        tmp_phase_bin);
                    
                    
                    tmp_is_valid_Q = ~arrayfun(@(x) isempty(x.median), tmp_flux_bin_stat);
                    tmp_phase_val = tmp_phase_val(tmp_is_valid_Q);
                    tmp_flux_med = [tmp_flux_bin_stat.median].';
                    tmp_flux_2575 = cat(1, tmp_flux_bin_stat.prctile_val);
                    tmp_flux_2575 = tmp_flux_2575(:, [7, 9]);
                    tmp_flux_err = abs(tmp_flux_2575 - tmp_flux_med);
                    
                    
                    errorbar(tmp_ax_hdl, tmp_phase_val, tmp_flux_med, ...
                        tmp_flux_err(:, 1), tmp_flux_err(:, 2));
                    tmp_ax_hdl.XLim = [-pi, pi];
                    tmp_ax_hdl.YLim = flux_lim;
                    grid(tmp_ax_hdl, 'on');
                    tmp_ax_hdl.Title.String = sprintf('r \\in [%.2f, %.2f] \\mum (#%d)', ...
                        r_bin_edge(i_bin), r_bin_edge(i_bin+1), numel(tmp_phase));
                    tmp_ax_hdl.XLabel.String = 'Phase (rad)';
                    tmp_ax_hdl.YLabel.String = 'Flux (RBC/s)';
                end
            end
            if ~dev_mode
                vis_folder = fullfile(tmp_vis_folder_r, 'flux_vs_phase');
                fig_fp = fullfile(vis_folder, 'ROI', sprintf('%s_%s', tmp_vis_folder_n, 'f_vs_pl_bin_by_r.png'));
                fun_print_image_in_several_formats(fig_hdl, fig_fp);
                delete(fig_hdl);
            end
        end
    end
end