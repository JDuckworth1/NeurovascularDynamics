

% function [is_near_midline_Q] = Thomas_generate_hemisphere_mask(vsl_mask,vsl_im,neuron_mask)
function [is_near_midline_Q,mix_gaussian_model] = Thomas_generate_hemisphere_mask(vsl_mask,vsl_im,neuron_mask)
%% Parameters
neuron_mask = neuron_mask.rim;
% Remove connected components with less than 50 pixels (you might want to
% change this number according to the images)
min_cc_size = 100;

% Assume that the midline of the brain is roughly along the x-axis. The two
% centers of the gaussian should be well-separated along the y direction. 
min_hemisphere_center_dist_y = 350; %Was 200

% Midline probability half-width
% Classify the pixels with classification ambiguity to be the pixels near
% the midline:
midline_prob_half_width = 0.25;

%%
mask_size = size(vsl_mask);
assert(all(mask_size == size(vsl_im)) & all(mask_size == size(neuron_mask)), 'Inconsistent mask size');

vsl_mask_clean = bwareaopen(vsl_mask, min_cc_size);
% imshowpair(vsl_mask_clean, vsl_mask);
% Get pixel coordinates
vsl_ind = find(vsl_mask_clean);
[vsl_sub_y, vsl_sub_x] = ind2sub(mask_size, vsl_ind);
vsl_xy = cat(2, vsl_sub_x, vsl_sub_y);
% Use 2 gaussians to fit the position of the vessel mask pixels 
options = statset('MaxIter',1000);
% options = statset('Display','final');
mix_gaussian_model = fitgmdist(vsl_xy, 2,'Options',options);
fprintf('Obtain the following mixture of gaussians statistics:\n');
disp(mix_gaussian_model);

if abs(diff(mix_gaussian_model.mu(:, 2))) < min_hemisphere_center_dist_y
    warning('The centers of gaussians is less than %d pixels along the y-direction', min_hemisphere_center_dist_y);
    is_near_midline_Q = NaN;
    return %Added
end
% Classification
prob_in_mixture = mix_gaussian_model.posterior(vsl_xy);
is_in_mixture_1 = prob_in_mixture(:, 1) > prob_in_mixture(:, 2);
% Generate mask for each hemisphere
mix_1_mask = false(mask_size);
mix_1_mask(vsl_ind(is_in_mixture_1)) = true;
mix_1_mask_cvxh = bwconvhull(mix_1_mask);

mix_2_mask = false(mask_size);
mix_2_mask(vsl_ind(~is_in_mixture_1)) = true;
mix_2_mask_cvxh = bwconvhull(mix_2_mask);
% Generate the mask for the midline: 
[im_y, im_x] = ndgrid(1 : mask_size(1), 1 : mask_size(2));
prob_1 = mix_gaussian_model.posterior(cat(2, im_x(:), im_y(:)));
is_near_midline_Q = abs(prob_1(:, 1) - 0.5) < midline_prob_half_width;
is_near_midline_Q = reshape(is_near_midline_Q, mask_size);
%% Visualization
fig_hdl = figure;
% fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,2];
ax_hdl = subplot(2,2,1);
imagesc(ax_hdl, vsl_mask);
% colormap('gray');
ax_hdl.DataAspectRatio = [1,1,1];
gmPDF = @(x,y) arrayfun(@(x0, y0) pdf(mix_gaussian_model, [x0, y0]), x, y);
hold(ax_hdl, 'on');
fc_hdl = fcontour(ax_hdl, gmPDF, [ax_hdl.XLim, ax_hdl.YLim]);
fc_hdl.LineWidth = 2;
ax_hdl.Title.String = 'Vessel mask with Gaussian mixture contour';

ax_hdl_2 = subplot(2,2,2);
imagesc(ax_hdl_2, vsl_im .* uint8(mix_1_mask_cvxh));
colormap(ax_hdl_2, 'gray');
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.Title.String = 'Masked image 1';

ax_hdl_3 = subplot(2,2,3);
imagesc(ax_hdl_3, vsl_im .* uint8(mix_2_mask_cvxh));
colormap(ax_hdl_3, 'gray');
ax_hdl_3.DataAspectRatio = [1,1,1];
ax_hdl_3.Title.String = 'Masked image 2';

ax_hdl_4 = subplot(2,2,4);
% imagesc(ax_hdl_4, vsl_im .* uint8(is_near_midline_Q));
imagesc(ax_hdl_4, vsl_im + uint8(is_near_midline_Q));
colormap(ax_hdl_4, 'gray');
ax_hdl_4.DataAspectRatio = [1,1,1];
ax_hdl_4.Title.String = 'Mask for midline';
end