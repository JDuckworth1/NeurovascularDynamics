function graph_str = fun_skeleton_to_graph_2D(voxel_list, mask_size)
% This function computes the graph from the skeleton. 
% Input:
%   voxel_list: N-by-1 double precision array of skeleton voxel indices
%   mask_size: 2-by-1 double precision array. size of the mask(skeleton
%   array)
% Output: graph_str with fields:
%   node:structure with fields:
%         pos_ind: linear index position of the node voxel in the mask 
%         num_voxel, num_cc(#node): each node can have more than 1 voxels
%         label: N-by-1 double precision array. specify the label of the
%         node voxel in pos_ind;
%         cc_ind: cell array, each contains the linear index of the voxel
%         in each node
%         num_link: number of link that this node attached
%         connected_link_label: labels of links that this node joins.
%         map_ind_2_label: N-by-1 sparse double precision matrix for
%         mapping the voxel linear index to its label.
%   num: structure with fields
%         mask_size(2-by-1), skeleton_voxel(number of skeleton voxel),
%         mask_size_pad( = mask_size +2, i.e. pad one layer of 0 on both
%         size), blk_vol(number of voxels in the block, or mask),
%         neighbor_add_pad(8-by-1 double precision array for finding the
%         26 neighbors, in the padded array), and block_voxel_pad.
%   link, isopoint and endpoint: strucutres with similar fileds with node. 
%
% Implemented by Xiang Ji on 03/25/2019
% Adapted from fun_skeleton_to_graph
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
num.neighbor_add_pad = sub2ind(num.mask_size_pad, tmp1(:), tmp2(:));
num.neighbor_add_pad = num.neighbor_add_pad - num.neighbor_add_pad(5);
num.neighbor_add_pad(5) = [];

num.neighbor_add = sub2ind(num.mask_size, tmp1(:), tmp2(:));
num.neighbor_add = num.neighbor_add - num.neighbor_add(5);
num.neighbor_add(5) = [];
%% Generate 1D sparse matrix representation of the skeleton array
[pos_1, pos_2] = ind2sub(num.mask_size, voxel_list);
num.skeleton_voxel = numel(voxel_list);
voxel_idx_padded = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1); 
sp_skl = sparse(voxel_idx_padded, ones(num.skeleton_voxel, 1), ...
    1:num.skeleton_voxel, num.block_voxel_pad, 1);
clear pos_1 pos_2
%% Classify voxels
% list of neighbor voxel index in voxel_list for each voxel. 
l_voxel_neighbor_idx = full(sp_skl(bsxfun(@plus, voxel_idx_padded', num.neighbor_add_pad)));
% Sort the neighbor for link tracking in the next section
l_voxel_neighbor_idx = sort(l_voxel_neighbor_idx, 1, 'descend');
l_num_neighbor = sum(l_voxel_neighbor_idx>0,1);
% clear sp_skl
% List of endpoint index in the voxel list
l_ep_idx = find(l_num_neighbor == 1);
num_endpoints = numel(l_ep_idx);
% List of node index in the voxel list
l_nd_idx = find(l_num_neighbor > 2);
% List of link point logical index in the voxel list
l_lp_Q = (l_num_neighbor == 2);

l_isop_Q = (l_num_neighbor == 0);
%% Track along links
l_neighbor_label_of_nd_idx = l_voxel_neighbor_idx(:, l_nd_idx);
l_neighbor_label_of_nd_idx = l_neighbor_label_of_nd_idx(l_neighbor_label_of_nd_idx>0);
l_neighbor_label_of_nd_idx = unique(l_neighbor_label_of_nd_idx);
% The following vector contains both the endpoint indices and link point
% indices. 
l_neighbor_link_label_of_nd_idx = setdiff(l_neighbor_label_of_nd_idx, l_nd_idx);
% The start tacking points are the union of the link points in the neighbor
% to the node points, or the end points ( in case that the links have two
% endpoints)
l_link_start_voxel_idx = union(l_neighbor_link_label_of_nd_idx, l_ep_idx);
t_num_link_start_voxel = numel(l_link_start_voxel_idx);
% Label the voxels need to visit
voxel_unvisited = true(num.skeleton_voxel,1);
voxel_unvisited(l_nd_idx) = false;
voxel_unvisited(l_isop_Q) = false;
num_unvisited_points = nnz(l_lp_Q) + num_endpoints;
t_link_idx = zeros(1000,1);
t_start_search_idx = 1;
t_num_cc = 0;
link_cc.PixelIdxList = cell(num.skeleton_voxel, 1);
% while num_unvisited_points > 0
t_start_point_idx = 0;
while t_start_point_idx < t_num_link_start_voxel
    % Find the starting voxel list index in the voxel list - for links
    for t_start_point_idx = t_start_search_idx : t_num_link_start_voxel
        t_num_cc_voxel = 0;   
        t_current_id = l_link_start_voxel_idx(t_start_point_idx);
        if voxel_unvisited(t_current_id)
            t_start_search_idx = t_start_point_idx + 1;
            keep_tracking = true;
            break;
        end
    end
    while keep_tracking
        keep_tracking = false;
        % Add the current voxel to the connected component voxel list 
        t_num_cc_voxel = t_num_cc_voxel + 1;
        % MATLAB can extend the length of the list automatically if
        % t_num_cc_voxel is larger than the initialized length. 
        t_link_idx(t_num_cc_voxel) = t_current_id;
        voxel_unvisited(t_current_id) = false;
        % Get the neighbors of the current voxel and pick the ONE hasn't
        % been visited. 
        t_neighbor_idx = l_voxel_neighbor_idx(:, t_current_id);
        for tmp_idx = 1 : 2
            tmp_id = t_neighbor_idx(tmp_idx);
            if tmp_id > 0
                if voxel_unvisited(tmp_id)
                    keep_tracking = true;
                    t_current_id = tmp_id;
                    break
                end
            else
                break;
            end
        
        end
    end
    if t_num_cc_voxel > 0
        t_num_cc = t_num_cc + 1;
        num_unvisited_points = num_unvisited_points - t_num_cc_voxel;
        link_cc.PixelIdxList{t_num_cc} = voxel_list(t_link_idx(1:t_num_cc_voxel));
    end
end
link_cc.PixelIdxList = link_cc.PixelIdxList(1:t_num_cc)';
link_cc.NumObjects = t_num_cc;
%% Construct graph
graph_str = struct;
graph_str.num = num;

% Link information
link_length_list = cellfun(@length, link_cc.PixelIdxList');
% Should be very carefully when trying to delete some links. Links,
% endpoints and nodes should be modified in the same time. 
graph_str.link.num_cc = link_cc.NumObjects;
graph_str.link.cc_ind = link_cc.PixelIdxList';
graph_str.link.pos_ind = cat(1, graph_str.link.cc_ind{:});
graph_str.link.num_voxel = numel(graph_str.link.pos_ind);
graph_str.link.num_voxel_per_cc = link_length_list;

% Transpose to make a column vector
graph_str.link.label = repelem(1:graph_str.link.num_cc, graph_str.link.num_voxel_per_cc)';
graph_str.link.map_ind_2_label = sparse(graph_str.link.pos_ind, ...
    ones(graph_str.link.num_voxel,1), ...
    graph_str.link.label, ...
    graph_str.num.block_voxel,1);

% Endpoint information (Notice that endpoints are a subset of link points)
graph_str.endpoint.pos_ind = voxel_list(l_ep_idx);
graph_str.endpoint.link_label = full(graph_str.link.map_ind_2_label(graph_str.endpoint.pos_ind));
graph_str.endpoint.num_voxel = num_endpoints;
graph_str.endpoint.map_ind_2_label = sparse(graph_str.endpoint.pos_ind, ones(graph_str.endpoint.num_voxel, 1),...
    1 : graph_str.endpoint.num_voxel, graph_str.num.block_voxel, 1);

% Isolated point (Just record for completeness...)
graph_str.isopoint.pos_ind = voxel_list(l_isop_Q);
graph_str.isopoint.num_voxel = nnz(graph_str.isopoint.pos_ind);
% Isolated loop
if num_unvisited_points > 0
    warning('Isolated loops exist in the skeleton. Ignored...');
    loop_cc = fun_cc_in_sparse_matrix(voxel_list(voxel_unvisited), mask_size);
    graph_str.isoloop.cc_ind = loop_cc.PixelIdxList';
    graph_str.isoloop.num_cc = loop_cc.NumObjects;
end

%% Correcting node information
tmp_mask = false(num.mask_size);
tmp_mask(voxel_list(l_nd_idx)) = true;
node_cc = bwconncomp(tmp_mask);

% Node information
graph_str.node.num_voxel = nnz(l_nd_idx);
graph_str.node.num_cc = node_cc.NumObjects;
% Transpose the cell array
graph_str.node.cc_ind = node_cc.PixelIdxList';
graph_str.node.pos_ind = cat(1, node_cc.PixelIdxList{:});
graph_str.node.num_voxel_per_cc = cellfun(@length, graph_str.node.cc_ind);

% Transpose the vector to make it a column vector
graph_str.node.label = repelem(1:graph_str.node.num_cc, graph_str.node.num_voxel_per_cc)';
graph_str.node.map_ind_2_label = sparse(graph_str.node.pos_ind, ...
    ones(graph_str.node.num_voxel,1), ...
    graph_str.node.label, ...
    graph_str.num.block_voxel,1);
%% Connect nodes with links
% Since the order of the voxel indices has been changed, the linear
% position of each voxel cannot be extracted by the logical array anymore.
% Need to re-calculate or find a list indices map. 
[pos_1, pos_2] = ind2sub(num.mask_size, graph_str.link.pos_ind);
link_pos_ind_pad = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1); 
[pos_1, pos_2] = ind2sub(num.mask_size, graph_str.node.pos_ind);
node_pos_ind_pad = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1); 
% Be careful about the construction of the map look up table!!!
link_map_ind_pad_2_lable = sparse(link_pos_ind_pad, ones(graph_str.link.num_voxel, 1), ...
    graph_str.link.label, num.block_voxel_pad,1);

graph_str.node.num_link = zeros(graph_str.node.num_cc,1);
graph_str.node.connected_link_label = cell(graph_str.node.num_cc,1);
graph_str.link.connected_node_label = zeros(graph_str.link.num_cc,2);
graph_str.link.num_node = zeros(graph_str.link.num_cc,1);
% for each voxel in the node list, look at the 26 neighbors in the padded
% link point look up table. 
for tmp_node_idx = 1 : graph_str.node.num_voxel
    tmp_node_label = graph_str.node.label(tmp_node_idx);

    node_neighbor_link_label = full(link_map_ind_pad_2_lable(node_pos_ind_pad(tmp_node_idx) + num.neighbor_add_pad));
    node_neighbor_link_label = node_neighbor_link_label(node_neighbor_link_label>0);
    node_neighbor_link_label = unique(node_neighbor_link_label);
    num_neighbor_seg = length(node_neighbor_link_label);
    for tmp_neighbor_idx = 1 : num_neighbor_seg
        tmp_l_voxel_neighbor_idx = node_neighbor_link_label(tmp_neighbor_idx);
        % Add information to node
        graph_str.node.num_link(tmp_node_label) = graph_str.node.num_link(tmp_node_label) + 1;
        graph_str.node.connected_link_label{tmp_node_label}(graph_str.node.num_link(tmp_node_label)) = tmp_l_voxel_neighbor_idx;
        % Add information to link        
        graph_str.link.num_node(tmp_l_voxel_neighbor_idx) = graph_str.link.num_node(tmp_l_voxel_neighbor_idx) + 1;
        graph_str.link.connected_node_label(tmp_l_voxel_neighbor_idx, graph_str.link.num_node(tmp_l_voxel_neighbor_idx)) = tmp_node_label;
    end
end

end
