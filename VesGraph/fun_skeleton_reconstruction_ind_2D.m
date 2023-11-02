function recon_voxel_ind = fun_skeleton_reconstruction_ind_2D(pos_ind, block_size, recon_strel)
% 
% Input: 
%   pos_ind: N-by-1 numerical array, indice of the voxel in the block 
%   block_size: 2-by-1 numerical array, the size of the block 
%   strel_array: 2 dimensional logical array. The structure element for
%   reconstruction. 
% Output:
%   recon_mask: 2 dimensional logical array. 
% 
% This function is equivalent to imdilate(skeleton, strel). This function
% is faster than imdilate if the structure element is large( say more than
% 5x5x5)

% Compute the position of the voxel in the structure element
strel_size = size(recon_strel);
% num_block_voxel = prod(block_size);
[pos1, pos2] = ind2sub(strel_size, find(recon_strel));
% Pad the array
strel_r = max((size(recon_strel) - 1)/2);
block_size_pad = block_size + strel_r * 2;
recon_strel_add_pad = sub2ind(block_size_pad, pos1 + strel_r, pos2 + strel_r);
recon_strel_add_pad = recon_strel_add_pad - recon_strel_add_pad(ceil(length(recon_strel_add_pad)/2));
% Subscript of the centerline voxel in the block 
[pos1, pos2] = ind2sub(block_size, pos_ind);
% Indices of the centerline voxel in the padded block
pos_ind_pad = sub2ind(block_size_pad, pos1 + strel_r, pos2 + strel_r);
% Indices of the reconstructed mask in the padded block 
recon_voxel_ind_pad = bsxfun(@plus, pos_ind_pad, recon_strel_add_pad');
% Subscripts of the reconstructed mask in the padded block 
[pos1, pos2] = ind2sub(block_size_pad, recon_voxel_ind_pad(:));
% Subscripts of the reconstructed mask in the bounding box of the original
% block in the padded block 
select_Q = (pos1 > strel_r & pos1 <= (block_size(1) + strel_r)) & ...
    (pos2 > strel_r & pos2 <= (block_size(2) + strel_r));
% Subscript of the reconstructed mask on 
pos1 = pos1 - strel_r;
pos2 = pos2 - strel_r;
pos1(~select_Q) = nan;
pos2(~select_Q) = nan;
recon_voxel_ind = sub2ind(block_size, pos1, pos2);
recon_voxel_ind = reshape(recon_voxel_ind, size(recon_voxel_ind_pad));
end
