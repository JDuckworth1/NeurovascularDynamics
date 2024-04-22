function est_bg_str = fun_DKLab_TBRL_background_estiamtion(vsl_im, est_sig_int_n, est_snr, vis_Q)

if nargin < 4
    vis_Q = false;
end

est_bg_str = struct;
est_bg_str.est_sig_int_n = est_sig_int_n;
est_bg_str.est_snr = est_snr;

vsl_im_mean = mean(vsl_im, 3);
im_size = size(vsl_im_mean);
% Mainly for removing the signal in the first row - from line scan? 
vsl_im_mean_v =  movmedian(mean(vsl_im_mean, 2), 5);
vsl_im_mean_v_n = vsl_im_mean_v ./ max(vsl_im_mean_v);

est_bg_str.avg_im = vsl_im_mean;
est_bg_str.avg_v_int = vsl_im_mean_v;
est_bg_str.avg_v_int_n = vsl_im_mean_v_n;

int_r = [];

while numel(int_r) < 2 && est_sig_int_n < 1
    est_bg_str.est_sig_int_n = est_sig_int_n;
    int_r = fun_find_zeros_by_linear_interpolation(1 : numel(vsl_im_mean_v_n), ...
        vsl_im_mean_v_n - est_sig_int_n);
    int_r = round(int_r);
    est_sig_int_n = est_sig_int_n * 1.1;
end

est_bg_str.bg_row_ind = int_r;
if isscalar(int_r)
    warning('Only one side of the average image has average pixel value lower than the esimated threshold');
    if int_r <= im_size(1)/2
        bg_sample_row_ind = 1: 1 : int_r;
    else
        bg_sample_row_ind = int_r : 1 : vsl_im_mean(1);
    end
elseif numel(int_r) >= 3
    warning('The average intensity profile has more than 2 points corssing the estimated threshold');
    dist_r0 = diff(int_r);
    [~, max_idx] = max(dist_r0);
    bg_sample_row_ind = [1 : int_r(max_idx), int_r(max_idx+1), im_size(1)];
    % Use the two with greater gradient
%     int_diff = [0; diff(vsl_im_mean_v_n)];
%     int_r_diff = int_diff(int_r);
else
    bg_sample_row_ind = [1:int_r(1), int_r(2):im_size(1)];   
end
est_bg_str.bg_row_ind = bg_sample_row_ind;
est_bg_str.sig_row_ind = setdiff(1 : im_size(1), bg_sample_row_ind, 'stable');

est_bg_str.bg_mean = mean(vsl_im_mean_v(bg_sample_row_ind));
est_bg_std = std(single(vsl_im(bg_sample_row_ind, :, :)), 0, 3);
est_bg_str.bg_std = double(mean(est_bg_std, 'all'));
est_bg_str.est_bg_int = est_bg_str.bg_mean + est_bg_str.bg_std * est_bg_str.est_snr;
est_bg_str.est_sig_avg_int = mean(est_bg_str.avg_v_int(est_bg_str.sig_row_ind));
if vis_Q
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
    plt_hdl = plot(ax_hdl, vsl_im_mean_v);
    hold(ax_hdl, 'on');
    line_hdl = line(ax_hdl, [1, im_size(1)], [est_bg_str.est_bg_int, est_bg_str.est_bg_int], ...
        'Color', 'r');
    legend(ax_hdl, 'Average vertical profile', ...
        sprintf('Estiamted background (SNR %.1f)', est_snr));    
end

end