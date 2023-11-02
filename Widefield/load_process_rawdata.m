% load_process_rawdata.m

%Code to load wide-field tiff images, extract df/f(x,t) and spectra
%Processes before-stim, stim, and after-stim time periods separately (deter
%= 'b', 'p', 'a').
%Written for separate tiff images acquired by Prime95 camera.
%Alternate image loading functions can be used for loading tiff stacks.
%Written to process one imaging session at a time (4-8 experiments). Masks
%can be created once and used for all experiments in the imaging session
%(brain does not move during the session).
%Uses adi reader from:
%https://github.com/JimHokanson/adinstruments_sdk_matlab to load LabChart
%adi files with behavior and camera frame information.

%Written by Thomas Broggini, Jacob Duckworth @UCSD
%Contact: Jacob jaduckwo@ucsd.edu

clc;clear;close all;
%Create files struct containing the list of files
%Change directories to the folder containing data
% files = dir('10.Oct.2023*\*meters.mat') %Get files by date (imaging session)
files = dir('*TB013120M4*\*meters.mat') %Get files by animal
deter = 'b' %options b=beforpuff p=puff a=afterpuff

for n=1:length(files)
    
    cd(files(n).folder) %% Set Parameter for Data Loading

    %NrOCh = '1' %options 1 = one channel, 2G = green of dual, 2R = red, 2 for Ratiometricmeasurements. of dual this is inversed when using SMART to control the LED exposure
    str = 'v' %Analyze vessels or neurons v/n [v]; %% Loading Preparations
    close all
    rmpath(genpath('C:\chronux_2_12')); %Get rid of CHRONUX in path so we can use the matlab findpeaks function.
    data_folder = files(n).folder;
    cd(data_folder);
    file   = dir('*_Parameters.mat')
    filename = extractBefore(files(n).name,'_Parameters.mat');
    fname = erase(file.name,'_Parameters.mat');
    str1=["dualCH","PwrG80R","RedGreen","PowerG80R","dualLED"];
    str2=["G80_","singleCH","_Single_","5Power","SingleCH","G05"];

    [NrOCh] = fun_get_NrOCh(str,fname,str1,str2); %Determine number of channels for loading

    load(files(n).name) %Load parameters file

    if ~isfile('Mask1.tif')
        [im_cell] = fun_load_ims(NrOCh,0,2500); %Load images for mask1

        im_data = cat(3, im_cell{:});
        clear im_cell
        diffusedImage = imdiffusefilt(im_data,'NumberOfIterations',5);
        target = std(double(diffusedImage),0,3);
        clear diffusedImage
        addpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis')
        contrast = adapthisteq(normfunc(double(target)));
        figure
        imagesc(contrast); daspect([1,1,1])
        cd(files(n).folder)
        imwrite(contrast, 'Mask1.tif');
        clear im_data
    end

    target = imread('Mask1.tif');
    if ~isfile('rim.mat')
        figure; imagesc(target);
        daspect([1,1,1]);
        rim = roipoly();
        close
        save('rim.mat','rim')
    end
    load('rim.mat')
    target=double(target).*rim;
    cd(data_folder);
    addpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\adi_reader\')
    disp(fname)

    sensordata_path = fullfile(['\\dk-server.dk.ucsd.edu\jaduckwo\RVLMSensorData\','\',uigetfile('\\dk-server.dk.ucsd.edu\jaduckwo\RVLMSensorData\*.adicht','Select an Sensor File')])
    Sensorprams=adi.readFile(sensordata_path)

    if contains(Sensorprams.channel_names{7},'LED');
        tmp_stimLED = Sensorprams.getChannelByName('Stim LED');
        SensorData.stimLED=tmp_stimLED.getData(1);
    end
    tmp_p95b=Sensorprams.getChannelByName('Prime95 Output');
    try
    tmp_leftp = Sensorprams.getChannelByName('left puffer');
    catch
        tmp_leftp = [];
    end
    tmp_rightp = Sensorprams.getChannelByName('right puffer');
    tmp_leftp = Sensorprams.getChannelByName('left puffer');
    tmp_x=Sensorprams.getChannelByName('X');
    tmp_y = Sensorprams.getChannelByName('Y');
    tmp_z = Sensorprams.getChannelByName('Z');
    SensorData.p95b=tmp_p95b.getData(1);

    SensorData.rightp=tmp_rightp.getData(1);
    SensorData.leftp = tmp_leftp.getData(1);
    SensorData.x=tmp_x.getData(1);
    SensorData.y=tmp_y.getData(1);
    SensorData.z=tmp_z.getData(1);
    % Parameters
    timeLimits = [0 Sensorprams.channel_specs(1,4).n_samples]; % seconds
    frequencyLimits = [0 20]; % Hz
    % Index into signal time region of interest
    tmp_p95b_ROI = SensorData.p95b(:);
    sampleTime = Sensorprams.channel_specs(1,4).dt; % seconds
    startTime = 0; % seconds
    minIdx = ceil(max((timeLimits(1)-startTime)/sampleTime,0))+1;
    maxIdx = floor(min((timeLimits(2)-startTime)/sampleTime,length(tmp_p95b_ROI)-1))+1;
    tmp_p95b_ROI = tmp_p95b_ROI(minIdx:maxIdx);
    [p,f]=pspectrum(tmp_p95b_ROI,seconds(sampleTime),'FrequencyLimits',frequencyLimits);
    figure('units','normalized','outerposition',[0 0 .5 1]);
    plot(f,p)
    
    [wind,~] = ginput(2);
    N = 200 ;
    fi = linspace(wind(1),wind(2)) ;
    yi = interp1(f,p,fi);

    %Get maximum
    [~,fpeak1] = max(yi) ;
    acq_frequency = fi(fpeak1)
    clear fpeak1 wind
    im_acq_freq = acq_frequency;
    if NrOCh ~= '1'
        acq_frequency=acq_frequency/2;
    end
    toplot.rate = acq_frequency;
    close
    if contains(Sensorprams.channel_names{6},'LED')
        figure
        plot(SensorData.stimLED)
    end
    aprams=struct('folder', data_folder, 'imfolder', im_folder, 'senorprams', Sensorprams, 'camdiff', sampleTime);

    %% Load Prime Data Prepeartions
    %Need to write a more generalized function that can accept any possible
    %LabChart adi fields. For now this works (needs to be updated depending
    %on experiment).

    close all
    rmpath(genpath('C:\chronux_2_12'));
    rmpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
    clearvars -except files transferfiles files_tranferred filename n str toplot wave rim pside recprams data_folder fname im_mask im_folder SensorData Sensorprams acq_frequency sampleTime fpeak1 aprams deter NrOCh im_acq_freq BW nROIs
    [~,idx] = findpeaks(diff(SensorData.p95b),'MinPeakHeight',1.5);
    if isfield(SensorData,'rightp')
        [~,pridx] = findpeaks(diff(SensorData.rightp),'MinPeakHeight',1);
    end
   
    if contains(fname,'Mod')
        [~,stimidx] = findpeaks(diff(SensorData.stimLED),'MinPeakHeight',0.009);
    else
        [~,stimidx] = findpeaks(diff(SensorData.stimLED),'MinPeakHeight',1);
%         [~,stimidx] = findpeaks(diff(SensorData.leftp),'MinPeakHeight',1);
    end

if isfield(SensorData,'rightp')
    aprams.rpuffdiff=pridx;
end
%Defines with images of the Trial to load
if deter == 'b'
    tmp_start = 1 ;%Always start with 1 like Matlab
    tmp_end = stimidx(1)-1;
elseif deter == 'p'
    tmp_start = stimidx(1);
    tmp_end = stimidx(length(stimidx))-1;
elseif deter == 'a'
    tmp_start = stimidx(length(stimidx));
    tmp_end = idx(length(idx));
end
sampleTime = Sensorprams.channel_specs(1,4).dt; % seconds
tmp_im_start=round(tmp_start*sampleTime*im_acq_freq)
tmp_im_end=round(tmp_end*sampleTime*im_acq_freq )-1
im_sec_list = tmp_im_start:tmp_im_end;
num_im = numel(im_sec_list);
im_cell = cell(num_im, 1);
if NrOCh == '2'
    im_cell2 = cell(num_im, 1);
end

toplot.fname = fname;
if deter == 'b';
    toplot.fname=[toplot.fname,'svd_beforepuff'];
elseif deter == 'p';
    toplot.fname=[toplot.fname,'svd_puff'];
elseif deter == 'a';
    toplot.fname=[toplot.fname,'svd_afterpuff'];
end
if NrOCh == '2'
    toplot.fname=[toplot.fname,'_ratio'];
end
if str == 'n'
    toplot.fname=[toplot.fname,'_neurons'];
end
cd(aprams.folder);

if str == 'v'
    save([toplot.fname,'.mat'])
elseif str =='n'
    save([toplot.fname,'.mat'])
end
disp('Done Saving')
%% Perform Loading loop

tic
cd(data_folder);
[im_cell] = fun_load_ims(NrOCh,tmp_im_start,tmp_im_end); %Load images for whole experiment time period

im_data = cat(3, im_cell{:});
clear im_cell
toc

%% Skeletonize Prime Data for Vessel Calcium
%Finalize mask with photoshop or gimp before this section.
if str == 'v'
    cd(data_folder)
    %Pre-process mask (correct for RGB)
    if ~isfile('Mask.tif')
        [mask_file,mask_file_path] = uigetfile('*.tif','Select an Mask File'),
        im_mask = ~logical(imread(fullfile(mask_file_path,mask_file)));
        if nnz(im_mask ) > nnz(~im_mask) == 1
            im_mask=~im_mask;
        end
        imwrite(im_mask, 'Mask.tif');
    else
        mask_file='Mask.tif';
        mask_file_path=data_folder;
        %im_mask = logical(imread(fullfile(mask_file_path,mask_file)));
        im_mask = imread(fullfile(mask_file_path,mask_file));
        if ndims(im_mask) ==3
            im_mask = im_mask(:,:,1)>10;
        else
            im_mask = logical(im_mask);
        end
        if nnz(im_mask ) > nnz(~im_mask) %Check if we need to invert the mask.
            im_mask=~im_mask;
        end
        imwrite(im_mask, 'Mask.tif');
    end

    if ~isfile('Mask.tif')
        diffusedImage = imdiffusefilt(im_data,'NumberOfIterations',5);
        target = std(double(diffusedImage),0,3);
        clear diffusedImage
        contrast = adapthisteq(normfunc(double(target)));
        imwrite(contrast, 'Mask1.tif');
    end
    % %Then use Photoshop to finalize the map Manualy
    toplot.target=mean(im_data,3);
    if size(im_mask,3)>1
        im_mask=squeeze(im_mask(:,:,1));
    end
    imwrite(im_mask, 'Mask.tif');

    %Ensure mask is properly oriented on the brain.
    figure
    A = toplot.target;

    %Optional
%     Nmax = 1000; %get Nmax biggest entries
%     [ Avec, Ind ] = sort(A(:),1,'descend');
%     A(Ind(1:Nmax)) = mean(A,'all');

    overlay = cat(3,normfunc(A').*255,0.75.*255.*im_mask');
    overlay(:,:,3) = overlay(:,:,2); overlay(:,:,2) = overlay(:,:,1); overlay(:,:,1) = 0;
    imshow(uint8(overlay));imshow(uint8(overlay));
    toplot.mask=im_mask;

    % Extract data after skeletonization
    % For the meaning of the following fields of the structure, refer the
    % document of fun_get_skeleton_neighbor_ind
    % Use the vessel graph functions (written by Xiang Ji)

    addpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\XianCode\XianCode')

    opt.downsample_rate = 2; %Downsample rate for skeleton pixels
    opt.imdilate_disk_r = 1; %Dilate skeleton pixels by this * local vessel radius
    opt.min_kept_cc_num_pixel = 35; % You would need a larger value for higher resolution image 35 for Prime wide-field imaging
    opt.rm_neighbor_out_of_mask_Q = true;
    opt.vis_mask_Q = true; % Show the cleaned up mask
    [skel_neighbor_ind, toplot.mask, vsl_graph] = fun_get_skeleton_neighbor_ind(im_mask, opt);
    imwrite(toplot.mask, 'Mask.tif');
    [skel_neighbor_ind, ~] = fun_get_skeleton_neighbor_ind(toplot.mask, opt);
    if str =='v'
        % ### wave is computed here it 2D frames x sceletal element###
        if NrOCh == '2'
            wave = fun_get_skeleton_neighbor_stat_from_image_stack(im_data, skel_neighbor_ind, 'median');
            wave2 = fun_get_skeleton_neighbor_stat_from_image_stack(im_data2, skel_neighbor_ind, 'median');
            wave = ((wave./mean(wave, 2))./(wave2./mean(wave2,2)))-1;
            clear wave2
        else
            wave = fun_get_skeleton_neighbor_stat_from_image_stack(im_data, skel_neighbor_ind, 'median');
            wave = bsxfun(@rdivide ,bsxfun(@minus, wave, mean(wave, 2)),mean(wave, 2));%calculate dF/F
        end
        % ###
        [skel_neighbor_ind_unique, tmp_unique_idx, ~] = unique(cat(2, skel_neighbor_ind{:}), 'stable');
        recon_mask = false(size(im_mask));
        recon_mask(skel_neighbor_ind_unique) = true;
        num_skel_ind_unique = numel(skel_neighbor_ind_unique);
        tmp_skel_label = repelem(1 : numel(skel_neighbor_ind), cellfun(@numel, skel_neighbor_ind));

        skel_label_unique = tmp_skel_label(tmp_unique_idx);
        toplot.skel_label = skel_label_unique;
        toplot.mask_ind = skel_neighbor_ind_unique;
        toplot.mask_size = size(im_mask);
        disp('Finish extracting data for skeleton pixels');
        [U,S,V]=svd(wave,0); %Space-time SVD
        for i=1:size(S,1)
            lambda(i)=S(i,i)^2;
        end
        figure
        plot(log(lambda));
        xlim([0 100]);
        toplot.fname = fname;
        if deter == 'b';
            toplot.fname=[toplot.fname,'svd_beforepuff'];
        elseif deter == 'p';
            toplot.fname=[toplot.fname,'svd_puff'];
        elseif deter == 'a';
            toplot.fname=[toplot.fname,'svd_afterpuff'];
        end
        if NrOCh == '2'
            toplot.fname=[toplot.fname,'_ratio'];
        end
        if isfile([toplot.fname,'_wave.h5']);
            delete([toplot.fname,'_wave.h5']);
        end
        h5create([toplot.fname,'_wave.h5'],'/wave',[size(wave,1), size(wave,2)],'Datatype','double');
        h5write([toplot.fname,'_wave.h5'],'/wave', wave, [1 1],[size(wave,1), size(wave,2)]);
        clear wave
    end
end

%% Preprocess Neuronal Data
if str == 'n'
%Use the mask to extract data from the image stack
cd(data_folder)
im_data_rs = reshape(im_data, numel(im_data(:, :, 1)), [])';
% %only do this for binniing
% %im_data = imresize3(im_data,[toplot.mask_size(1)/2 toplot.mask_size(2)/2 size(im_data,3)],'Method','nearest');
toplot.target = mean(im_data,3);
if ~isfile('Mask.tif')
    [mask_file,mask_file_path] = uigetfile('*.tif','Select an Mask File')
    im_mask = ~logical(imread(fullfile(mask_file_path,mask_file)));
    if nnz(im_mask ) > nnz(~im_mask) == 1;
        im_mask=~im_mask;
    end
    imwrite(im_mask, 'Mask.tif');
else
    mask_file='Mask.tif';
    mask_file_path=data_folder;
    %             im_mask = logical(imread(fullfile(mask_file_path,mask_file)));
    im_mask = imread(fullfile(mask_file_path,mask_file));
    if ndims(im_mask) ==3
        im_mask = im_mask(:,:,1)>10;
    else
        im_mask = logical(im_mask);
    end
    if nnz(im_mask ) > nnz(~im_mask) == 1
        im_mask=~im_mask;
    end
    imwrite(im_mask, 'Mask.tif');
end
conn = bwconncomp(im_mask);
if conn.NumObjects>20
    opt.downsample_rate = 2;
    opt.imdilate_disk_r = 1;
    opt.min_kept_cc_num_pixel = 35; % You would need a larger value for higher resolution image 35 for Prime
    opt.rm_neighbor_out_of_mask_Q = true;
    opt.vis_mask_Q = true; % Show the cleaned up mask
    [~, im_mask] = fun_get_skeleton_neighbor_ind(im_mask, opt);
end
toplot.mask = im_mask;
imwrite(toplot.mask, 'Mask.tif');

%Cut the outside black to reduce the pixelnumber
if exist('rim.mat','file') == 1 %was == 2
    load('rim.mat')
end
if exist('rim','var')
    im_mask=(~toplot.mask.*rim);
else
    rim = roipoly();
    im_mask=(~toplot.mask.*rim);
end
figure
A = toplot.target;
Nmax = 1000; % get Nmax biggest entries
[ Avec, Ind ] = sort(A(:),1,'descend');
max_values = Avec(1:Nmax);
[ind_row, ind_col] = ind2sub(size(A),Ind(1:Nmax)); % fetch indices
A(ind_row, ind_col)=mean(A,'all');
addpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis')
overlay = cat(3,normfunc(A).*255,0.75.*255.*im_mask);
overlay(:,:,3) = overlay(:,:,2); overlay(:,:,2) = overlay(:,:,1); overlay(:,:,1) = 0;
imshow(uint8(overlay));imshow(uint8(overlay)); 
daspect([1,1,1]);
% end
%Does it fit?
overlay = cat(3,normfunc(toplot.target').*255,0.75.*255.*im_mask');
overlay(:,:,3) = overlay(:,:,2); overlay(:,:,2) = overlay(:,:,1); overlay(:,:,1) = 0;
figure;
imshow(uint8(overlay));
%Extract Pixels and generate a Matrix
im_mask_ind = find(imresize(im_mask,1));
wave = double(im_data_rs(:, im_mask_ind))';
toplot.mask_ind = im_mask_ind;
num_skel_ind_unique = numel(toplot.mask_ind);
skel_label_unique = 1 : num_skel_ind_unique;
toplot.skel_label = skel_label_unique;
toplot.mask_size = size(im_mask);
%wave = double(skel_data_ori);
%clear im_data_rs
clear im_data
wave = bsxfun(@rdivide ,bsxfun(@minus, wave, mean(wave, 2)),mean(wave, 2)); %calculate dF/F
[num_pixel,num_frame] = size(wave);
[U,S,V]=svd(wave,0);
for i=1:size(S,2)
    lambda(i)=S(i,i)^2;
end
figure
plot(log(lambda));
xlim([0 100]);
toplot.fname = fname;
if deter == 'b';
toplot.fname=[toplot.fname,'svd_beforepuff'];
elseif deter == 'p';
toplot.fname=[toplot.fname,'svd_puff'];
elseif deter == 'a';
toplot.fname=[toplot.fname,'svd_afterpuff'];
end
if str == 'n'
toplot.fname=[toplot.fname,'_neurons'];
end
if isfile([toplot.fname,'_wave.h5']);
    delete([toplot.fname,'_wave.h5']);
end

cd(data_folder)
h5create([toplot.fname,'_wave.h5'],'/wave',[size(wave,1), size(wave,2)],'Datatype','double');
h5write([toplot.fname,'_wave.h5'],'/wave', wave, [1 1],[size(wave,1), size(wave,2)]);
clear wave
disp('Done Saving')
end
%% Parameters
if str=='v'

% Update half-bandwidth according to the rounded p value
if length(V) < length(U)
    sig_modes = round(length(V)/8);
else
    sig_modes = length(U);
end
Un=single(U(:,1:sig_modes));
Sn=single(S(1:sig_modes,1:sig_modes));
Vn=single(V(:,1:sig_modes));
if str == 'n';
    clear U S V
end

%Get stim frequency from file name
ext=extractBetween(fname,'_0.','Hz');
if isempty(ext)
    ext=extractBetween(fname,'_0,','Hz');
end
if isempty(ext)
    prompt = 'What is the stimfreqency? ';
    fpeak2 = input(prompt)
    fpeak3=fpeak2*2;
    toplot.f_peak(2)=fpeak2;
    toplot.f_peak(3)=fpeak3;
else
    fpeak2=str2double(erase(['0.',ext{1}],'_'));
    if size(ext,1)>2;
        fpeak3=str2double(erase(['0.',ext{2}],'_'));
        fpeak4=str2double(erase(['0.',ext{3}],'_'));
        fpeak5=fpeak2*2;
        toplot.f_peak(4)=fpeak4;
        toplot.f_peak(5)=fpeak5;
    elseif size(ext,1)>1;
        fpeak4=fpeak2*2;
        fpeak3=str2double(erase(['0.',ext{2}],'_'));
        toplot.f_peak(4)=fpeak4;
    else
        fpeak3=fpeak2*2;
    end
toplot.f_peak(2)=fpeak2;
toplot.f_peak(3)=fpeak3;
end
toplot.Delta_f = 0.02; % Hz. denoted as half-bandwidth 

num_pixel = size(Un,1);
num_frame= size(Vn,1);
padding_ratio = 2; 
num_frame_pad = (2 ^ ceil(log2(num_frame))) * padding_ratio;
toplot.pad = num_frame_pad;
p = round(num_frame * toplot.Delta_f / toplot.rate); % (num_frame / acq_frequency = acquisition time, i.e. T in spectral methods)
toplot.Delta_f = p * toplot.rate / num_frame;
disp(['Bandwidth = ', num2str(toplot.Delta_f), ' Hz'])
toplot.num_tapers = 2 * p - 1; %Max number of tapers without distortion
[slep,~] = dpss(num_frame, p, toplot.num_tapers);

%The puff toplot.f_peak(1) is the beforepuff vasomotion freq
if deter == 'p' || deter == 'a'
    btoplot = load((append(data_folder,'\',filename,'svd_beforepuff.mat')),'toplot'); 
    btoplot = btoplot.toplot;
    toplot.f_peak(1) = btoplot.f_peak(1);
end
%% Perform Spectral FFT
addpath(genpath('C:\chronux_2_12'))

ntapers = [(toplot.num_tapers+1)/2,toplot.num_tapers];
Fs = toplot.rate;
nfft = toplot.pad;

tapers = dpsschk(ntapers,size(Vn,1),Fs);
[f,findx] = getfgrid(Fs,nfft,[0,Fs/2]);

%FFT
taperedFFT = complex(zeros(length(f),ntapers(2),size(Vn,2)));
for i = 1:size(Vn,2) %multi-taper FT of temporal space-time svd components
    J = mtfftc(Vn(:,i),tapers,nfft,Fs);
    J = J(findx,:,:);
    taperedFFT(:,:,i) = J;
end

scores = Un*Sn; %How to reconstruct wave 
clear Un Sn i J tapers
rmpath(genpath('C:\chronux_2_12'));

%% Look at Spectral Power & Map the S_tot back to the mask

S_tot = zeros(size(scores,1),size(taperedFFT,1),'single');
if str=='v'
    parfor k = 1:ntapers(2)
        z = scores*squeeze(taperedFFT(:,k,:))';
        S_tot = S_tot + conj(z).*z;
    end
else
    for k = 1:ntapers(2)
        tic
        n = size(scores,1);

        blksize = 22000;
        gtaperedFFT=gpuArray(single(complex(squeeze(taperedFFT(:,k,:))')));

        for i=0:(n/blksize-1)
            rows = i*blksize + (1:blksize);
            gscores=gpuArray(scores(rows,:));
            S_tot(rows, :) = S_tot(rows, :) + gather(conj(gscores*gtaperedFFT).*(gscores*gtaperedFFT));
        end

        if rows(end)
            rows = (rows(end)+1):n;
            gscores=gpuArray(scores(rows,:));
            gtaperedFFT=gpuArray(single(complex(squeeze(taperedFFT(:,k,:))')));
            S_tot(rows, :) = S_tot(rows, :) + gather(conj(gscores*gtaperedFFT).*(gscores*gtaperedFFT));
            
        end

        toc
    end
    clear('gscores','gtaperedFFT');
end
S_tot = S_tot./ntapers(2);
clear Fs ntapers nfft tapers findx
[num_pixel] = size(toplot.mask_ind,1);
num_frame_pad=size(S_tot,2);
toplot.mpowr=mean(S_tot,1);
toplot.f = f;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(f, log(mean(S_tot, 1)), 'ko');
if deter == 'b'% && str == 'n'
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    hold on
    if str == 'v'
        for i = 1 : size(S_tot,1)
            plot(f, log(S_tot(i,:)));
        end
    else
        for i = 1 :100: size(S_tot,2)
            plot(f, log(S_tot(i,:)));
        end
    end
    plot(f, log(mean(S_tot, 1)), 'ko');
    subplot(2,1,2)
    toplot.mpowr=mean(S_tot,1);
    plot(log10(toplot.mpowr),'ko')
    %title(fname, 'Interpreter', 'none');
    prompt='Define Start and end of the peak search: ';
    [toplot.wind,~] = ginput(2);
    toplot.wind = fix(toplot.wind);

    toplot.pwr = [];
    hold off
    %We transform back
    tic
    if str == 'v'
        S_tot = S_tot(toplot.skel_label,:);
    end
    [~, toplot.pwr] = max(S_tot(:, toplot.wind(1):toplot.wind(2)), [], 2);
    toplot.pwr = toplot.f(toplot.pwr + toplot.wind(1) - 1);
    toc
    disp(' Done');
end
    %Find the max using the delta function
    if str == 'v' && deter == 'b'
        [~, fpeak1] = max(toplot.mpowr(:, toplot.wind(1):toplot.wind(2)), [], 2);
        fpeak1 = toplot.f(fpeak1 + toplot.wind(1) - 1);

        toplot.f_peak=[];
        toplot.f_peak(1)=fpeak1;
        k=extractBetween(fname,'_0.','Hz');
        fpeak2=str2double(erase(['0.',k{1}],'_'));
        if size(k,1)>1
            fpeak3=str2double(erase(['0.',k{2}],'_'));
        else
            fpeak3=fpeak2*2;
        end
        toplot.f_peak(2)=fpeak2;
        toplot.f_peak(3)=fpeak3;
    end
end

%% Save Data
toplot.startframe = im_sec_list(1);
toplot.endframe = im_sec_list(end);
toplot.fname = fname;
if deter == 'b';
toplot.fname=[toplot.fname,'svd_beforepuff'];
elseif deter == 'p';
toplot.fname=[toplot.fname,'svd_puff'];
elseif deter == 'a';
toplot.fname=[toplot.fname,'svd_afterpuff'];
end
if NrOCh == '2'
    toplot.fname=[toplot.fname,'_ratio'];
end
if str == 'n'
toplot.fname=[toplot.fname,'_neurons'];
end
%cd(aprams.folder);
% cd(append('\\dk-server.dk.ucsd.edu\jaduckwo\GalProject_Analysis\',fname));
if str == 'v'
save([toplot.fname,'.mat'], 'toplot' , 'recprams' ,'Sensorprams','SensorData', 'aprams', 'Vn','scores');
% save([toplot.fname,'.mat'], 'toplot' , 'recprams' ,'Sensorprams','SensorData', 'aprams');
elseif str =='n'
%     save([toplot.fname,'.mat'], 'toplot' , 'recprams' , 'Sensorprams','SensorData', 'aprams' ,'Vn', 'rim'); %'-append');
    save([toplot.fname,'.mat'], 'toplot' , 'recprams' , 'Sensorprams','SensorData', 'aprams' , 'rim'); %'-append');
end
disp('Done Saving')
close all
%save([toplot.fname,'.mat'], 'toplot' , 'recprams' , 'Sensorprams','SensorData', 'aprams', 'rim');
end