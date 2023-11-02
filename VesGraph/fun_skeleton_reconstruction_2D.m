function mask = fun_skeleton_reconstruction_2D(skeleton_voxel, skeleton_radius, block_size)

% Classify skeleton voxel
skeleton_radius = round(skeleton_radius);
radius_list = unique(skeleton_radius);
% Binary reconstruction
mask = false(block_size);
for tmp_idx = 1 : numel(radius_list)
    tmp_radius = radius_list(tmp_idx);
    if tmp_radius >= 1
        tmp_ind = fun_skeleton_reconstruction_ind_2D(skeleton_voxel(skeleton_radius == tmp_radius),...
            block_size, strel('disk', double(tmp_radius)).Neighborhood);
        mask(tmp_ind) = true;
    end
end
    
end

