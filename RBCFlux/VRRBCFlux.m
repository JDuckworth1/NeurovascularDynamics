classdef VRRBCFlux
    
    properties
        
    end
    
    
    methods
        function obj = VRRBCFlux()
            
        end
    end
%% Measurement
    methods(Static)
        function im = concatenate_frames(im)
            [~, ~, num_frame] = size(im);
            if num_frame > 1
                im_cell = cell(num_frame, 1);
                for i = 1 : num_frame
                    im_cell{i} = im(:, :, i);
                end
                im = cat(2, im_cell{:});
            end
        end
        
        function [vsl_radius_um, varargout] = line_radius_estimation(int_profile, ...
                est_sig_min_int, est_bg_mean, pxl_size_um, vis_Q)
            if nargin < 5
                vis_Q = false;
            end
            % Initialization
            varargout = {};
            int_profile_str = struct;
            int_profile_str.r_zeros = nan(1, 2);
            int_profile_str.radius_um = nan;
            int_profile_str.threshold = nan;
            
            %
            int_profile_str.int = int_profile;
            int_profile_str.x_pxl = (1 : numel(int_profile));
            
            try
                seg_local_max_int = findpeaks(int_profile, int_profile_str.x_pxl, 'MinPeakHeight', est_sig_min_int, ...
                    'MinPeakDistance', 5);
                seg_threshold_int = (median(seg_local_max_int) - est_bg_mean) / 2 + est_bg_mean;
            catch
                seg_threshold_int = nan;
            end
            
            % Find the zero-crossing position
            seg_zero_pxl = [];
            while numel(seg_zero_pxl) < 2 && seg_threshold_int < intmax('uint16')
                seg_zero_pxl = fun_find_zeros_by_linear_interpolation(int_profile_str.x_pxl,...
                    int_profile - seg_threshold_int);
                seg_threshold_int = seg_threshold_int * 1.1;
            end
            
            int_profile_str.threshold = seg_threshold_int;
            if numel(seg_zero_pxl) >= 2
                vsl_radius_um = (seg_zero_pxl(end) - seg_zero_pxl(1))/2 * pxl_size_um;
                int_profile_str.r_zeros_pxl = seg_zero_pxl([1, end]);
            else
                vsl_radius_um = nan;
                int_profile_str.r_zeros_pxl = [nan, nan];
                warning('The intensity profile should cross zero at least twice');
            end
            
            if nargout > 1
                int_profile_str.radius_um = vsl_radius_um;
                int_profile_str.r_zeros = int_profile_str.r_zeros_pxl .* pxl_size_um;
                varargout{1} = int_profile_str;
            end
            
            % Visualization
            if vis_Q
                fig_hdl = figure;
                ax_hdl = axes(fig_hdl);
                plot(ax_hdl, int_profile_str.x_pxl, int_profile);
                hold(ax_hdl, 'on');
                line(ax_hdl, ax_hdl.XLim, [seg_threshold_int, seg_threshold_int], ...
                    'Color', 'g');
            end
        end
        
        function [fig_hdl] = vis_image_section(im_r, im_f, vis_sec, pad_sec)
            if nargin < 4
                pad_sec = 50;
            end
            [~, num_lines] = size(im_r);
            vis_sec_list = max(vis_sec - pad_sec, 1) : min(vis_sec + pad_sec, num_lines);
            vis_im_r = im_r(:, vis_sec_list);
            vis_im_f = im_f(:, vis_sec_list);
            
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1.5, 1];
            ax_hdl_1 = subplot(2,1,1);
            ax_hdl_2 = subplot(2,1,2);
            imagesc(ax_hdl_1, vis_im_r);
            imagesc(ax_hdl_2, vis_im_f);
            ax_hdl_2.XTickLabel = arrayfun(@(x) num2str(x, '%d'), ax_hdl_2.XTick + vis_sec_list(1) - 1, ...
                'UniformOutput', false);
            ax_hdl_1.XAxis.Visible = 'off';
        end
        
%         function [r_trace, str_cell] = estimate_radius_for_each_line(im_r, ...
%                 est_sig_int, est_bg_mean, pxl_size_um)
%             % im_r: double 2D array, each column is a line scan image
%             [~, num_lines] = size(im_r);
%             r_trace = nan(num_lines, 1);
%             str_cell = cell(num_lines, 1);
%             tic_r = tic;
%             for i = 1 : num_lines
%                 tmp_int = im_r(:, i);
%                 [r_trace(i), str_cell{i}] = VRRBCFlux.line_radius_estimation(...
%                     tmp_int, est_sig_int, est_bg_mean, pxl_size_um);
%                 if mod(i, 5000) == 0
%                     fprintf('Finish processing %.2f%% of the lines. Elpased time is %.2f seconds\n', ...
%                         (i/num_lines) * 100, toc(tic_r));
%                 end
%             end
%         end
        
        function vsl_mask = compute_vessel_mask(vsl_im, est_th, min_cc_size, disk_r, sph_r)
            vsl_im_stack_size = size(vsl_im);
            vsl_mask = false(vsl_im_stack_size);
            t_tic = tic;
            for i_frame = 1 : vsl_im_stack_size(3)
                tmp_im = vsl_im(:, :, i_frame);
                tmp_mask = imbinarize(tmp_im, 'global');
                if min(tmp_im(tmp_mask)) < est_th
                    tmp_mask = tmp_im > tmp_th;
                end
                tmp_mask = bwareaopen(tmp_mask, min_cc_size);
                tmp_mask = imclose(tmp_mask, strel('disk', disk_r));
                % Keep the largest connected component
%                 tmp_cc = bwconncomp(tmp_mask);
%                 tmp_cc_sz = cellfun(@numel, tmp_cc.PixelIdxList);
%                 [~, max_idx] = max(tmp_cc_sz);
%                 tmp_mask = false(size(tmp_mask));
%                 tmp_mask(tmp_cc.PixelIdxList{max_idx}) = true;                
                vsl_mask(:, :, i_frame) = tmp_mask;
                if mod(i_frame, 1000) == 0
                    fprintf('Finish processing %d%%\n', ...
                        round(100 * i_frame / vsl_im_stack_size(3)));
                end
            end
            vsl_mask = imclose(vsl_mask, strel('sphere', sph_r));
            fprintf('Finish computing vessel mask. Elapsed time is %.2f seconds\n', ...
                toc(t_tic));
        end
        
        function result_cell = compute_rbc_cc(rbc_im, vsl_mask, est_bg_fraction, ...
                min_rbc_size)
            
            if ~isfloat(rbc_im)
                rbc_im = single(rbc_im);
            end
            vis_Q = false;
            
            num_frame = size(rbc_im, 3);
            assert(all(size(rbc_im) == size(vsl_mask)));
            result_cell = cell(num_frame, 1);
            t_tic = tic;
             for i_frame = 1 : num_frame
            %for i_frame = 10 : 15
                result_cell{i_frame} = VRRBCFlux.single_frame_RBC_segmentation(...
                    rbc_im(:, :, i_frame), vsl_mask(:, :, i_frame), ...
                    est_bg_fraction, min_rbc_size, vis_Q);
                if mod(i_frame, 1000) == 0
                    fprintf('Finish processing %d%%\n', ...
                        round(100 * i_frame / num_frame));
                end
            end
            fprintf('Finish computing RBC connected components. Elapsed time is %.2f seconds\n', ...
                toc(t_tic));
        end
        
        function result_str = single_frame_RBC_segmentation(rbc_im, vsl_mask, ...
                est_bg_fraction, min_rbc_size, vis_Q)
            if nargin < 5
                vis_Q = false;
            end
            
            if ~isfloat(rbc_im)
                rbc_im = single(rbc_im);
            end
            result_str = struct;
            % Estimate background
            int_in_vsl_mask = rbc_im(vsl_mask);
            num_pxl_in_mask = numel(int_in_vsl_mask);
            result_str.num_pxl_in_mask = num_pxl_in_mask;
            int_sort = sort(int_in_vsl_mask, 'ascend');
            result_str.masked_int_med = int_sort(round(num_pxl_in_mask / 2));
            result_str.masked_bg_std = std(int_sort(1 ...
                : round(num_pxl_in_mask * est_bg_fraction)));
            est_th = result_str.masked_int_med + 3 * result_str.masked_bg_std;
            
            rbc_mask = vsl_mask & (rbc_im > est_th);
            rbc_cc = bwconncomp(rbc_mask);
            cc_sz = cellfun(@numel, rbc_cc.PixelIdxList);
            cc_valid_Q = cc_sz > min_rbc_size;
            rbc_mask(cat(1, rbc_cc.PixelIdxList{~cc_valid_Q})) = false;
            
            result_str.est_th = est_th;
            result_str.num_cc = nnz(cc_valid_Q);
            result_str.cc_ind = rbc_cc.PixelIdxList(cc_valid_Q);
            result_str.rbc_mask = rbc_mask;
            
            if vis_Q
                rbc_cc.PixelIdxList = rbc_cc.PixelIdxList(cc_valid_Q);
                rbc_cc.NumObjects = result_str.num_cc;
                tmp_rbc_lm = labelmatrix(rbc_cc);
                fig_hdl = figure;
                ax_hdl_1 = subplot(2,1,1);
                imagesc(ax_hdl_1, rbc_im);
                ax_hdl_2 = subplot(2,1,2);
                imagesc(ax_hdl_2, tmp_rbc_lm);
            end
        end
        
        function result_str = analyze_single_frame_rbc_cc(rbc_seg_cell, ...
                rbc_im)
            t_tic = tic;
            if ~isfloat(rbc_im)
                rbc_im = single(rbc_im);
            end
            im_size = size(rbc_im, [1,2]);
            num_frame = size(rbc_im, 3);
            assert(num_frame == numel(rbc_seg_cell));
            rbc_im_c = VRRBCFlux.concatenate_frames(rbc_im);
            rbc_mask = cellfun(@(x) x.rbc_mask, rbc_seg_cell, 'UniformOutput', false);
            rbc_mask = cat(3, rbc_mask{:});
            rbc_mask_c = VRRBCFlux.concatenate_frames(rbc_mask);
            % Recompute connected components
            rbc_cc = bwconncomp(rbc_mask_c);
            
            result_str.image_size = rbc_cc.ImageSize;
            result_str.num_cc = rbc_cc.NumObjects;
            result_str.cc_ind = rbc_cc.PixelIdxList;
            result_str.cc_int = cellfun(@(x) rbc_im_c(x), rbc_cc.PixelIdxList, ...
                'UniformOutput', false);
            result_str.num_pxl = cellfun(@numel, result_str.cc_ind).';
            [result_str.int_mean, result_str.int_std, ...
                result_str.t_mean, result_str.t_std, ...
                result_str.y_mean, result_str.y_std] = ...
                deal(zeros(result_str.num_cc, 1));
            
            result_str.rbc_trace = zeros(result_str.image_size(2), 1);
            
            for i_cc = 1 : result_str.num_cc
                tmp_int = result_str.cc_int{i_cc};
                result_str.int_mean(i_cc) = mean(tmp_int);
                result_str.int_std(i_cc) = std(tmp_int);
                tmp_sub = fun_ind2sub(rbc_cc.ImageSize, result_str.cc_ind{i_cc});
                tmp_sub_mean = mean(tmp_sub, 1);
                tmp_sub_std = std(tmp_sub, 1, 1);
                result_str.y_mean(i_cc) = tmp_sub_mean(1);
                result_str.y_std(i_cc) = tmp_sub_std(1);
                result_str.t_mean(i_cc) = tmp_sub_mean(2);
                result_str.t_std(i_cc) = tmp_sub_std(2);
                
                tmp_t = max(1, round(result_str.t_mean(i_cc)));
                result_str.rbc_trace(tmp_t) = result_str.rbc_trace(tmp_t) + 1;
            end   
            fprintf('Finish analyzing RBC connected components. Elapsed time is %.2f seconds\n', ...
                toc(t_tic));
        end
        
        function [vsl_radius_um, int_str] = frame_radius_estimation(im_r, ...
                im_mask, min_snr, pxl_size_um)
            if ~isfloat(im_r)
                im_r = single(im_r);
            end
            
            bg_mask = imerode(~im_mask, strel('disk', 5));
            bg_pxl = im_r(bg_mask);
            bg_mean = mean(bg_pxl);
            bg_std = std(bg_pxl);
            est_sig_int = bg_mean + min_snr * bg_std;
            int_prof = median(im_r, 2);
            [vsl_radius_um, int_str] = VRRBCFlux.line_radius_estimation(...
                int_prof, est_sig_int, bg_mean, pxl_size_um);
        end
        
        function r_str = estimate_radius_for_each_frame(im_r, vsl_mask, ...
               min_snr, pxl_size_um)
            t_tic = tic;
            num_frame = size(im_r, 3);
            r_str = cell(num_frame, 1);
            for i = 1 : num_frame
                [~, r_str{i}] = VRRBCFlux.frame_radius_estimation(im_r(:, :, i), ...
                    vsl_mask(:, :, i), min_snr, pxl_size_um);
                if mod(i, 1000) == 0
                    fprintf('Finish processing %d%%\n', round(100 * i / num_frame));
                end
            end          
            r_str = cat(1, r_str{:});
            fprintf('Finish estimating radius for each frame. Elapsed time is %.2f seconds\n', ...
                toc(t_tic));
        end
        
        function [r_trace, varargout] = estimate_radius_for_each_line(im_r, vsl_mask, ...
               min_snr, pxl_size_um)
            t_tic = tic;
            num_frame = size(im_r, 3);
            [r_trace, r_est_str] = deal(cell(num_frame, 1));
            for i = 1 : num_frame
                [r_trace{i}, r_est_str{i}] = VRRBCFlux.frame_radius_estimation_per_line(im_r(:, :, i), ...
                    vsl_mask(:, :, i), min_snr, pxl_size_um);
                if mod(i, 1000) == 0
                    fprintf('Finish processing %d%%\n', round(100 * i / num_frame));
                end
            end          
            r_trace = cat(1, r_trace{:});
            if nargout > 1
                varargout{1} = cat(1, r_est_str{:});
            end
            fprintf('Finish estimating radius for each frame. Elapsed time is %.2f seconds\n', ...
                toc(t_tic));
        end
        
        function [r_um, r_est_str] = frame_radius_estimation_per_line(im_r, im_mask, ...
                min_snr, pxl_size_um)
            if ~isfloat(im_r)
                im_r = single(im_r);
            end
            % Background estimation
            bg_mask = imerode(~im_mask, strel('disk', 5));
            bg_pxl = im_r(bg_mask);
            bg_mean = mean(bg_pxl);
            bg_std = std(bg_pxl);
            est_sig_int = bg_mean + min_snr * bg_std;
            num_lines = size(im_r, 2);
            r_um = nan(num_lines, 1, 'single');
            r_est_str = cell(num_lines, 1);
            for i = 1 : num_lines
                [r_um(i), r_est_str{i}]= VRRBCFlux.line_radius_estimation(...
                    im_r(:, i), est_sig_int, bg_mean, pxl_size_um);
            end
            r_est_str = cat(1, r_est_str{:});
        end
    end
%% Analysis
    methods(Static)

    end
%% Visualization
    methods(Static)
        function [fig_hdl] = vis_r_f_traces(r_trace, flux_trace)
            
            
        end
    end
end