
clear; clc; close all;
%Find 
files =  dir('17.Oct.2023*\*_puff.mat') %For deter == 'p'
% files =  dir('17.Oct.2023*\*_beforepuff.mat') %For deter == 'b'
str = 'v';
deter = 'p';

%% Perform Loading
for file =1:length(files)

    clearvars -except files file notprocessed str deter
    cd(files(file).folder)

    %Load toplot mat
    data_folder = files(file).folder;
    cd(data_folder);
    filename = extractBefore(files(file).name,'_puff.mat');
    fname = erase(files(file).name,'_puff.mat');
    if deter == 'p'
        toplottmp = load(append(data_folder,'\',fname,'_puff.mat'));
    elseif deter == 'b'
        toplottmp = load(append(data_folder,'\',fname));
    end
    toplot = toplottmp.toplot;
    if deter == 'p' || deter == 'a'
        btoplot = load((append(data_folder,'\',filename,'_beforepuff.mat')),'toplot');
        btoplot = btoplot.toplot;
        toplot.f_peak(1) = btoplot.f_peak(1);
    end

    %% Set parameters

    im_mask=toplot.mask;
    tmp_mode = 1; %Plot first (dominant) space-freq SVD mode
    map = zeros(toplot.mask_size);
    im_size = size(toplot.mask);
    im_phase = zeros(size(im_mask));
    im_mag=zeros(size(im_mask));
    tmp_max_frm = 10;
    tmp_min_frm = -10;
    tmp_max_phase = pi/2; % Note bounds and if any data fall outside of bounds
    tmp_min_phase = -pi/2; % Note bounds and if any data fall outside of bounds
    tmp_max_mag = 0.025;
    tmp_min_mag = 0;
    tmp_max_pwr = log10(0.001);
    tmp_min_pwr = log10(0.00001);
    tmp_max_f0 = 0.2;
    tmp_min_f0 = 0.0;
    toplot.prams.max_frm = tmp_max_frm;
    toplot.prams.min_frm = -tmp_min_frm;
    toplot.prams.max_phase = tmp_max_phase;
    toplot.prams.min_phase = tmp_min_phase;
    toplot.prams.max_mag = tmp_max_mag;
    toplot.prams.min_mag = tmp_min_mag;
    toplot.prams.max_pwr = tmp_max_pwr;
    toplot.prams.min_pwr = tmp_min_pwr;
    toplot.prams.binsize = linspace(-pi,pi,360);
    ext=extractBetween(files(file).name,'0.','Hz');
    
    %% _____________ Phase

    tmp_peak_idx = 1; %This is the vasomotor frequency toplot.f_peak(1)
    toplot.phase=squeeze(angle(toplot.U(tmp_peak_idx, :, tmp_mode)));
    toplot.tphase=angle(sum(toplot.U(tmp_peak_idx, :, tmp_mode)));
    toplot.ophase = squeeze(toplot.phase-toplot.tphase);
    toplot.ophase = mod(toplot.ophase+pi,2*pi)-pi;

    map = zeros(toplot.mask_size);
    map(toplot.mask_ind) = toplot.ophase; %Map phase to spatial location of vessels
    rescaled_pwr = nan(toplot.mask_size);
    tmp_pixel_value = map(toplot.mask_ind);
    tmp_pixel_value = min(tmp_max_phase, max(tmp_min_phase, tmp_pixel_value) );
    rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

    cmap = colormap('jet'); %Use cyclic colormap for phase
    int_to_cmap = linspace(tmp_min_phase, tmp_max_phase, size(cmap,1));
    non_nan_ind = find(map);
    num_nonnan = numel(non_nan_ind);
    rbg_pwr_list = zeros(num_nonnan, 3);
    rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    rgb_pwr = ones(3, prod(im_size)) * 0.5; %This sets grayscale background
    rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
    rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
    toplot.rgb_phase = permute(rgb_pwr, [3, 2, 1]);
    clear i j pixcolor cmap int_to_cmap
    figure;
    imagesc(flip(toplot.rgb_phase,2));
    axis image;
    box off;
    colormap jet;
    caxis([tmp_min_phase tmp_max_phase]);
    daspect([1,1,1]);
    axis off

    %% _____________ Magnitude
    tmp_peak_idx = 1; %This is the vasomotor frequency toplot.f_peak(1)
    tmp_mg_data = squeeze(abs(toplot.U(tmp_peak_idx, :, tmp_mode)));
    map = zeros(toplot.mask_size);
    map(toplot.mask_ind) = tmp_mg_data;
    rescaled_pwr = nan(toplot.mask_size);
    tmp_pixel_value = map(toplot.mask_ind);

    tmp_pixel_value = min(tmp_max_mag, max(tmp_min_mag, tmp_pixel_value) );
    rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;

    cmap = colormap('jet');
    int_to_cmap = linspace(tmp_min_mag, tmp_max_mag, size(cmap,1));
    non_nan_ind = find(map);
    num_nonnan = numel(non_nan_ind);
    rbg_pwr_list = zeros(num_nonnan, 3);
    rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
    rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));
    rgb_pwr = ones(3, prod(im_size)) * 0.5; %This sets grayscale background
    rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
    rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
    toplot.rgb_mag = permute(rgb_pwr, [3, 2, 1]);
    clear i j pixcolor cmap int_to_cmap
    figure;
    imagesc(flip(toplot.rgb_mag,2));
    axis image;
    box off;
    colormap jet;
    caxis([tmp_min_mag tmp_max_mag]);
    daspect([1,1,1]);
    axis off

end
