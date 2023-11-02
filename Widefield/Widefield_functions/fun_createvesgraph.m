% fun_createvesgraph.m

function [toplot,vsl_graph,im_mask] = fun_createvesgraph(im_mask,toplot)

opt.downsample_rate = 2;
opt.imdilate_disk_r = 1;
opt.min_kept_cc_num_pixel = 35;
opt.rm_neighbor_out_of_mask_Q = true;
opt.vis_mask_Q = true; % Show the cleaned up mask
[~, im_mask, vsl_graph] = fun_get_skeleton_neighbor_ind(im_mask, opt);
[skel_neighbor_ind, ~] = fun_get_skeleton_neighbor_ind(toplot.mask, opt);
[skel_neighbor_ind_unique, tmp_unique_idx, ~] = unique(cat(2, skel_neighbor_ind{:}), 'stable');
tmp_skel_label = repelem(1 : numel(skel_neighbor_ind), cellfun(@numel, skel_neighbor_ind));
skel_label_unique = tmp_skel_label(tmp_unique_idx);
toplot.skel_label = skel_label_unique;
toplot.mask_ind = skel_neighbor_ind_unique;

end