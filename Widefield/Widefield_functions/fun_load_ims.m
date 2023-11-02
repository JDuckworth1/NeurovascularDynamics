% fun_load_mask1_ims.m

function [im_cell] = fun_load_ims(NrOCh,tmp_im_start,tmp_im_end)


im_sec_list = tmp_im_start:tmp_im_end;
num_im = numel(im_sec_list);
im_cell = cell(num_im, 1);
im_folder=dir('**\img_000001001_Default_000.tif');
im_folder = im_folder.folder

tic
warning('off', 'all');
for iter_im = 1 : num_im
    tmp_im_id = im_sec_list(iter_im);
    if mod(iter_im, 100) == 0
        fprintf(sprintf('img_channel000_position000_time%09d_z000.tif\n',tmp_im_id));%for micro-manager
    end

    tmp_fn = fullfile(im_folder, sprintf('img_%09d_Default_000.tif',tmp_im_id)); %for micro-manager
    %             tmp_fn = im_folder(iter_im).name;
    %             cd(im_folder(1).folder);
    if isfile(tmp_fn)
        if NrOCh == '1'
            %single color imaging
            im_cell{iter_im} = imread(tmp_fn);
        elseif NrOCh == '2R'
            %dual color imaging
            if mod(tmp_im_id,2) == 1 %1=green 0=red
                im_cell{iter_im} = imread(tmp_fn);
            end
        elseif NrOCh == '2G'
            if mod(tmp_im_id,2) == 0 %1=green 0=red
                im_cell{iter_im} = imread(tmp_fn);
            end
        elseif NrOCh == '2'
            if mod(tmp_im_id,2) == 1 %1=green 0=red
                im_cell{iter_im} = imread(tmp_fn);
            end
            if mod(tmp_im_id,2) == 0 %1=green 0=red
                im_cell2{iter_im} = imread(tmp_fn);
            end
        end
    else
        fprintf('The file does not exist. %s\n', tmp_fn);
    end
end
warning('on','all');

end


