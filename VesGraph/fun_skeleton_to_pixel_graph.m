function graph_str = fun_skeleton_to_pixel_graph(voxel_list, mask_size)
% This function construct the pixel graph from the skeleton. 
% Input:
%   voxel_list: N-by-1 double precision array of skeleton voxel indices
%   mask_size: 2-by-1 double precision array. size of the mask(skeleton
%   array)
% Output: 
%   graph_str: undirected graph strucutre generated by graph, where each
%   edge is a connection between skeleton pixel to its directly connected
%   neigbhoring skeleton pixel 
% Implemented by Xiang Ji on 12/13/2019

%% Initialization
if ~isvector(voxel_list)
    % If the input is the skeleton (2D logical array), convert it to voxel
    % list. 
    if nargin < 2
        mask_size = size(voxel_list);
    end
    voxel_list = find(voxel_list);
end

num.mask_size = mask_size;
num.mask_size_pad = num.mask_size + 2;
num.block_voxel = prod(num.mask_size);
num.block_voxel_pad = prod(num.mask_size_pad);
% 8 neighbors relative indices position array:
tmp1 = [1 2 3 1 2 3 1 2 3 ];
tmp2 = [1 1 1 2 2 2 3 3 3 ];
neighbor_idx_2_dist = [sqrt(2), 1, sqrt(2), 1, 1, sqrt(2), 1, sqrt(2)];
num.neighbor_add_pad = sub2ind(num.mask_size_pad, tmp1(:), tmp2(:));
num.neighbor_add_pad = num.neighbor_add_pad - num.neighbor_add_pad(5);
num.neighbor_add_pad(5) = [];

num.neighbor_add = sub2ind(num.mask_size, tmp1(:), tmp2(:));
num.neighbor_add = num.neighbor_add - num.neighbor_add(5);
num.neighbor_add(5) = [];
% Generate sparse matrix
[pos_1, pos_2] = ind2sub(num.mask_size, voxel_list);
num.skeleton_voxel = numel(voxel_list);
voxel_idx_padded = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1); 
sp_skl = sparse(voxel_idx_padded, ones(num.skeleton_voxel, 1), ...
    1:num.skeleton_voxel, num.block_voxel_pad, 1);
clear pos_1 pos_2
%% Construct the graph
% list of neighbor voxel index in voxel_list for each voxel. 
l_voxel_neighbor_idx = full(sp_skl(bsxfun(@plus, voxel_idx_padded', num.neighbor_add_pad)));
% map the neighbors to distance
[neighbor_idx, ~]= find(l_voxel_neighbor_idx);
dist_to_neighbor = neighbor_idx_2_dist(neighbor_idx);
graph_end_idx = l_voxel_neighbor_idx(l_voxel_neighbor_idx > 0);
l_num_neighbor = sum(l_voxel_neighbor_idx>0,1);
graph_start_idx = repelem(1 : num.skeleton_voxel, 1, l_num_neighbor)';
assert(numel(graph_end_idx) == numel(graph_start_idx), 'The number of starting points is different from the number of end points');
graph_end_ind = voxel_list(graph_end_idx);
graph_start_ind = voxel_list(graph_start_idx);

graph_str = graph(graph_start_ind, graph_end_ind, dist_to_neighbor);
end