function neighbor_ind_list = fun_skeleton_reconstruction_2D_cell(skl_ind, skeleton_radius, block_size)

% Classify skeleton voxel
skeleton_radius = round(skeleton_radius);
radius_list = unique(skeleton_radius);
% Binary reconstruction
neighbor_ind_list = cell(numel(skl_ind), 1);
for tmp_idx = 1 : numel(radius_list)
    tmp_radius = radius_list(tmp_idx);
    tmp_Q = (skeleton_radius == tmp_radius);
    if tmp_radius >= 1
        tmp_ind = fun_skeleton_reconstruction_ind_2D(skl_ind(tmp_Q),...
            block_size, strel('disk', double(tmp_radius)).Neighborhood);
        neighbor_ind_list(tmp_Q) = num2cell(tmp_ind, 2);
    end
end
% Remove NaN from the ind list for each pixel 
for iter_ind = 1 : numel(skl_ind)
    tmp_ind = neighbor_ind_list{iter_ind};
    is_nan_Q = isnan(tmp_ind);
    if any(is_nan_Q)
        tmp_ind = tmp_ind(~is_nan_Q);
        neighbor_ind_list{iter_ind} = tmp_ind;
    end
end
end

