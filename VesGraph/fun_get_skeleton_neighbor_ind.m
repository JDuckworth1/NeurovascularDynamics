function [skl_neighbor_ind, im_mask_sm, varargout] = fun_get_skeleton_neighbor_ind(im_mask, opt)
% fun_get_skeleton_neighbor_ind computes the neighbor pixel linear indices
% of the skeleton pixel of the input mask. 
% Input: 
%   im_mask : 2D logical array, mask of the image (vessel). If not logical
%   array, assume nonzero pixels to be true anc convert to logical array
%   automatically. 
%   opt : structure with fields: 
%       imdilate_disk_r : numerical scalar, integer. radius of the disk
%       structure element for smoothing the im_mask by imdilate. 
%       min_kept_cc_num_pixel : numerical scalar, integer, minimum size of
%       connected component used for extracting the skeleton. Parameter for
%       bwareaopen 
%       downsample_rate : numerical scalar, interger, skip every
%       downsample_rate pixels along the skeleton link segments in the
%       reconstruction. 
%       vis_mask_Q : logical scalar, if true, show both the input image
%       mask and the smoothed image mask
%       rm_neighbor_out_of_mask_Q : if true, remove all the pixels that are
%       in the neighbor of the skeleton but not in the original mask
% Output: 
%   skl_neighbor_ind : N-by-1 cell array. N is the number of skeleton
%   pixels. Each cell is a numerical vector of the neighbor pixels of each
%   sampled skeleton pixel. 
% 
% Implemented by Xiang Ji on April 1, 2019
skl_neighbor_ind = [];
if isempty(im_mask)
    return;
end
% Convert the mask to logical array if it is not
if ~islogical(im_mask)
    im_mask = logical(im_mask);
end
if ~any(im_mask)
    return;
end
%% Setting parameters
if nargin < 2
    % Parameters for smooth and clean the mask 
    opt = struct;
    opt.imdilate_disk_r = 0;
    opt.min_kept_cc_num_pixel = 0;
    opt.downsample_rate = 1;
    opt.rm_neighbor_out_of_mask_Q = true;
    opt.vis_mask_Q = false;
end
if ~isfield(opt, 'imdilate_disk_r')
    opt.imdilate_disk_r = 0;
end
if ~isfield(opt, 'min_kept_cc_num_pixel')
    opt.min_kept_cc_num_pixel = 0;
end
if ~isfield(opt, 'downsample_rate')
    opt.downsample_rate = 1;
end
if ~isfield(opt, 'rm_neighbor_out_of_mask_Q')
    opt.rm_neighbor_out_of_mask_Q = true;
end
if ~isfield(opt, 'vis_mask_Q')
    opt.vis_mask_Q = false;
end
%% Processing
if opt.imdilate_disk_r > 0
    im_mask_sm = imclose(im_mask, strel('disk', opt.imdilate_disk_r));
end
if opt.min_kept_cc_num_pixel > 0
    % Remove small connected components
    im_mask_sm = bwareaopen(im_mask_sm, opt.min_kept_cc_num_pixel);
end
if opt.vis_mask_Q
    figure;
    subplot(1,2,1)
    imagesc(im_mask);
    title('Input image mask');
    daspect([1,1,1]);
    subplot(1,2,2)
    imagesc(im_mask_sm);
    title('Smoothed image mask');
    daspect([1,1,1]);
end
% Convert the mask to skeleton ( centerline ) 
skel = bwskel(im_mask_sm);
% Use distance transform to estimate the radius of the vessel
im_dt = bwdist(~im_mask_sm);
% Convert skeleton to graph - for downsampling along the skeleton segments
vessel_graph = fun_skeleton_to_graph_2D(skel);
% Use half of the points in the link for computation
if opt.downsample_rate > 1
    used_skl_ind = cell(vessel_graph.link.num_cc, 1);
    % Round the downsample rate
    downsample_rate = ceil(opt.downsample_rate);
    
    for iter_link = 1 : vessel_graph.link.num_cc
        tmp_ind = vessel_graph.link.cc_ind{iter_link};
        used_skl_ind{iter_link} = tmp_ind(1:downsample_rate:end);
    end
    skel_ind = cat(1, used_skl_ind{:});
else
    skel_ind = vessel_graph.link.pos_ind;
end
% Get the radius of the skeleton pixel
skel_r = im_dt(skel_ind);
% Get the neigobhor pixel linear indices of each skeleton pixel 
im_size = size(im_mask_sm);
skl_neighbor_ind = fun_skeleton_reconstruction_2D_cell(skel_ind, skel_r, im_size);

num_skl = numel(skl_neighbor_ind);
for iter_pixel = 1 : num_skl
    tmp_ind = skl_neighbor_ind{iter_pixel};
    % Delete nan from the skeleton neighbor pixel list
    tmp_delete_Q = isnan(tmp_ind);
    if opt.rm_neighbor_out_of_mask_Q
       tmp_is_out_of_mask_Q = ~im_mask(tmp_ind); 
       tmp_delete_Q = tmp_delete_Q | tmp_is_out_of_mask_Q;
    end
    skl_neighbor_ind{iter_pixel} = tmp_ind(~ tmp_delete_Q);
end
is_emptyQ = cellfun(@isempty, skl_neighbor_ind);
skl_neighbor_ind = skl_neighbor_ind(~is_emptyQ);

if nargout > 1
    varargout{1} = vessel_graph;
end
end