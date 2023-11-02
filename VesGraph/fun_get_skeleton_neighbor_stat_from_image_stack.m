function skel_data = fun_get_skeleton_neighbor_stat_from_image_stack(im_data, skel_neighbor_ind, method)


num_im = size(im_data, 3);
num_pixel = numel(im_data(:, :, 1));
im_data = reshape(im_data, num_pixel, num_im);
im_data = im_data'; % Transpose to get num_im x num_pixel
num_skel_pixel = numel(skel_neighbor_ind);
skel_data = zeros(num_im, num_skel_pixel);
for iter_skel_pixel = 1 : num_skel_pixel %Iterate over segments
    tmp_ind = skel_neighbor_ind{iter_skel_pixel};
    % Extract data for all the neighbor pixels of the skeleton pixel. 
    % size(tmp_data) = num_neighbor_pixel x num_im
    tmp_data = double(im_data(:, tmp_ind))';
    switch method
        case 'mean'
            skel_data(:, iter_skel_pixel) = mean(tmp_data, 1, 'omitnan');
        case 'median'
            skel_data(:, iter_skel_pixel) = median(tmp_data, 1, 'omitnan');
        case 'max'
            skel_data(:, iter_skel_pixel) = max(tmp_data, [], 1, 'omitnan');
        case 'min'
            skel_data(:, iter_skel_pixel) = min(tmp_data, [], 1, 'omitnan');
        case 'prctile75'
            skel_data(:, iter_skel_pixel) = prctile(tmp_data, 75, 1);
        case 'random'
            skel_data(:, iter_skel_pixel) = datasample(tmp_data, 1, 1);
        case 'sum'
            skel_data(:, iter_skel_pixel) = sum(tmp_data, 1)'; %Sum over pixels in segment
        otherwise
            error('Unsupported method');
    end
end
skel_data = skel_data';
end