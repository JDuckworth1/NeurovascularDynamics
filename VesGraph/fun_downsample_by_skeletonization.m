function data_str = fun_downsample_by_skeletonization(im_data, im_mask, opt)

opt.sample_method = 'median';

data_str = struct;
disp('Debug');

data_str.skel_neighbor_ind = fun_get_skeleton_neighbor_ind(im_mask, opt);
% ### wave is computed here ###
data_str.sampled_data = fun_get_skeleton_neighbor_stat_from_image_stack(im_data, data_str.skel_neighbor_ind, opt.sample_method);
% ###   
[skel_ind_unique, tmp_unique_idx, ~] = unique(cat(2, data_str.skel_neighbor_ind{:}), 'stable');

recon_mask = false(size(im_mask));
recon_mask(skel_ind_unique) = true;
data_str.recon_mask = recon_mask;
data_str.num_recon_pixel = numel(skel_ind_unique);
data_str.num_sampled_skeleton_ind = numel(data_str.skel_neighbor_ind);
tmp_skel_label = repelem(1 : data_str.num_sampled_skeleton_ind,...
    cellfun(@numel, data_str.skel_neighbor_ind));
% skel_label_unique specify the way to copy the centerline value to the its
% neighbor pixels. For example, if skel_label_unique(a:b) = i, it means the
% a-th to b-th pixel at position skel_ind_unique(a:b) share the value of
% the i-th skeleton pixel
data_str.skel_ind_label = tmp_skel_label(tmp_unique_idx);



end