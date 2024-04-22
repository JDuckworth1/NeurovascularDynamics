clc;clear;close all;
DataManager = FileManager;

rui_folder_root = '';
rui_data_folder_list = {'JD230627M997 101123 flux measurement'};
%%
for raw_data_idx = 1

    %% Source folder
    raw_data_folder = rui_data_folder_list{raw_data_idx};
    folder_labels = strsplit(raw_data_folder);
    stack = strjoin(folder_labels(1:2), '_');
    data_info = struct;
    data_info.dataset_name = 'DKLab';
    data_info.stack = stack;
    data_info.source_folder = fullfile(rui_folder_root, raw_data_folder);
    %% Acquisition parameters
    data_info.parameters.linescan_rate_Hz = 7910 * 2; %approx equal to line period (us), although may vary between trials.
    data_info.parameters.pixel_size_um = 0.5;
    %%
    assert(isfolder(data_info.source_folder), 'The source folder does not exist');
    
    file_info = dir(data_info.source_folder);
    file_info = {file_info.name};
    [match_start, match_end] = regexp(file_info, "Line\w*.tif", 'once');
    is_matched_Q = cellfun(@(x) ~isempty(x) && x == 1, match_start);
    file_info = file_info(is_matched_Q);
    file_basename = cellfun(@(x) strsplit(x(5:end-4), '_'), file_info, 'UniformOutput', false);
    file_basename = cat(1, file_basename{:});
    file_basename = cellfun(@str2num, file_basename);
    [file_basename, s_idx] = sortrows(file_basename, 'ascend');
   
    data_info.roi_list = file_basename(:, 1);
    data_info.record_idx_list = file_basename(:, 2);
    data_info.source_filenames = file_info(s_idx).';
    
    data_info.num_files = numel(data_info.source_filenames);
    
    data_info.save_folder = DataManager.fp_raw_data_folder(data_info.dataset_name, data_info.stack);
    data_info.filename = sprintf('%s_%s_raw_data_info.mat', data_info.dataset_name, data_info.stack);
    data_info.filepath = fullfile(data_info.save_folder, data_info.filename);
    
    DataManager.write_data(data_info.filepath, data_info);
    %% Seperate radius and RBC flux measurements
    [data_info.plasma_filename, data_info.RBC_filename] = deal(cell(data_info.num_files, 1));
    data_info.data_size = nan(3, data_info.num_files);
    warning('off')
    for iter_file = 1 : data_info.num_files
        roi_idx = data_info.roi_list(iter_file);
        record_idx = data_info.record_idx_list(iter_file);
        
        tmp_tic_0 = tic;
        file_path = fullfile(data_info.source_folder, data_info.source_filenames{iter_file});
        save_folder = DataManager.fp_raw_data_folder(data_info.dataset_name, stack);
        
        im = DataManager.load_single_tiff(file_path);
        % Convert to uint16 - truncate at 0
        im = uint16(im);
        im_size = size(im);
        
        im_vessel = im(:, :, 2:2:end); %All vessel Cy5.5 images
        im_rbc = im(:, :, 1:2:end); %All CSFE RBC images         
               
        save_fn_base = sprintf('%s_%s_roi_%03d_%05d', data_info.dataset_name, ...
            stack, roi_idx, record_idx);
        % Radius data
        vsl_data = struct;
        vsl_data.dataset_name = data_info.dataset_name;
        vsl_data.stack = stack;
        vsl_data.roi_idx = roi_idx;
        vsl_data.record_idx = record_idx;
        vsl_data.im = im_vessel;
        vsl_data.fn = sprintf('%s_vsl.mat', save_fn_base);
        vsl_data.fp = fullfile(save_folder, vsl_data.fn);
        tmp_tic = tic;
        DataManager.write_data(vsl_data.fp, vsl_data);
        fprintf('Finish writing the vessel image stack, elapsed time is %.2f seconds.\n', ...
            toc(tmp_tic));
        
        data_info.image_filename{iter_file} = vsl_data.fn;
        data_info.data_size(:, iter_file) = im_size;
        %% RBC data
        ls_data = struct;
        ls_data.dataset_name = data_info.dataset_name;
        ls_data.stack = stack;
        ls_data.roi_idx = roi_idx;
        ls_data.record_idx = record_idx;
        ls_data.im = im_rbc;
        ls_data.fn = sprintf('%s_rbc.mat', save_fn_base);
        ls_data.fp = fullfile(save_folder, ls_data.fn);
        tmp_tic = tic;
        DataManager.write_data(ls_data.fp, ls_data);
        
        data_info.linescan_filename{iter_file} = ls_data.fn;
        
        fprintf('Finish writing the RBC flux data, elapsed time is %.2f seconds.\n', ...
            toc(tmp_tic));
        fprintf('Finish preprocessing ROI %d. Elapsed time is %.2f secodns.\n', ...
            roi_idx, toc(tmp_tic_0));
    end
    DataManager.write_data(data_info.filepath, data_info);
end