%ExtractPADiam.m
%Use to load 2p tif images, extract diameter using threshold in Radon
%space, and save results in struct for each depth.

%Load images from tif stacks, use TiRS to fit perimeter to penetrating vessel cross section images
%Note: user needs to edit code to process stim and rest trials
%% Load and organize image files.
clear; close all;
trialnum = 1; %This is the experiment session/day
%Base folder varies
base_folder = 'Y:\Data backup';
trials = readmatrix('Trials.xlsx','OutputType','char');
data_folder = [base_folder,'\',cell2mat(trials(trialnum))]

animalnum = extractAfter(cell2mat(trials(trialnum)),'WT_');
animal = [extractBefore(cell2mat(trials(trialnum)),' WT'),erase(animalnum,'_')] %For file naming. this includes experiment day and animal number

cd(data_folder);
files = dir('*001.tif');
files([files.bytes] > 1e10) = []; %Delete large files (not image sequence)
todel = zeros(length(files),1);
for i=1:length(files)
    if contains(files(i).name,'ROI') || contains(files(i).name,'test')
        todel(i) = true;
    end
end
files(logical(todel)) = [];

%% Process individual imaging runs
file = 1; 

[PA,depth1,depth2,~,pix_um] = fun_get_PAMag_depth(files(file).name);
rate = 7.25; %Hz, constant for all trials.

%Extract number of images
file_name = string([files(file).folder,'\',files(file).name])
file_name_py = strrep(file_name,' ','**'); %Work around for spaces in file names
outvars = pyrunfile(sprintf("TiffDimDetector.py %s",file_name_py),'numims');
ims = str2double(outvars.char)

%Use TiffImReader.py
%Can separate loading of depths to reduce memory requirements.
tic
outvarsShal = pyrunfile(sprintf("TiffImReader.py %s",file_name_py), 'imsout',r1 = int16(0), r2 = int16(ims), r3 = int16(2)); tmp1 = double(outvarsShal);
outvarsDeep = pyrunfile(sprintf("TiffImReader.py %s",file_name_py), 'imsout', r1 = int16(1), r2 = int16(ims + 1), r3 = int16(2)); tmp2 = double(outvarsDeep);
avg_ves1 = shiftdim(tmp1,1);
avg_ves2 = shiftdim(tmp2,1);
clearvars -except avg_ves1 avg_ves2 files file data_folder im_size loadims animal PA depth1 depth2 pix_um rate stim_str ims pix_um
toc

%% Visualize raw data

% figure
% for i=1:size(avg_ves1,3)
% %     imagesc(medfilt2(avg_ves2(:,:,i),[medFiltRange,medFiltRange]));
%     imagesc(avg_ves1(:,:,i));
%     daspect([1,1,1]);
%     str = sprintf('%.2f',i/7.25);
%     title(str);
%     pause()
% end

%% Fit perimeters to the data
for current_cell = 1:2
clearvars -except avg_ves1 avg_ves2 files file data_folder im_size loadims animal PA depth1 depth2 pix_um rate stim_str current_cell rate pix_um

if current_cell == 1
    avg_im = mean(avg_ves1,3);
    maskedImage_t = zeros(size(avg_ves1));
elseif current_cell == 2
    avg_im = mean(avg_ves2,3);
    maskedImage_t = zeros(size(avg_ves2));
end

%Define centroid for vessel of interest. This is used for defining lateral
%movement of vessel with depth (dx)
figure; imagesc(avg_im); daspect([1,1,1]); title('Drag ROI');
roi = drawrectangle;
rect = round(roi.Position);
figure; imagesc(avg_im); daspect([1,1,1]); title('Click Vessel Center');
[windx,windy] = ginput(1);

%% Threshold radon
if current_cell == 1
    imgSeqin = avg_ves1;
elseif current_cell == 2
    imgSeqin = avg_ves2;
end

[dEq,CircVec,imgSeq,vesselPerim,time] = PAThresholdRadon(0.4,0.2,imgSeqin,rect,rate,current_cell,pix_um);

figure
plot(time,dEq); title('Extracted equivalent diameter (um)')
% figure
% plot(time,CircVec);

% a = avg_MajAxis/2;
% b = avg_MinAxis/2;
% th_A = pi*a*b;
% h = ((a-b)^2)/((a+b)^2);
% th_P = pi*(a+b)*(1+3*h/(10+sqrt(4-3*h))); %Wiki approximation
% th_Circ = (4*th_A*pi)/(th_P^2) %Circularity for an ellipse of this shape
% ex_Circ = mean(CircVec)

%% Plot TiRS results
figure;
for i=2900:3000
    tmp_im = imgSeq(:,:,i);

    tmp_im = rescale(tmp_im);
    perimInds = find(vesselPerim(:,:,i));
    [prow,pcol] = ind2sub(size(tmp_im),perimInds);
    perimIndsFull = sub2ind(size(imgSeq(:,:,i)),prow,pcol);
    imagesc(tmp_im); daspect([1,1,1]); clim([0 1]);
    
    hold on
    plot(pcol,prow,'r','LineStyle','none','Marker','.','MarkerSize',10);
    hold off
    title(sprintf('Frame %.0f',i));
    pause();
end 

%% Save time series diameter fits
cd(data_folder); %Choose saving location.

%%% Necessary for files with rest and stim in same image stack.
% seriesLength = length(dEq);
% stim_diam = diam(1:seriesLength/2);
% rest_diam = diam((seriesLength/2 + 1):end);
% stim_dEq = dEq(1:seriesLength/2);
% rest_dEq = dEq((seriesLength/2 + 1):end);
%%%

allstats = struct();
allstats(1).name = files(file).name;
if current_cell == 1
    allstats(1).depth = depth1;
elseif current_cell == 2
    allstats(1).depth = depth2;
end
allstats(1).folder = files(file).folder;
allstats(1).RadondEq = dEq;
allstats(1).time = time;
allstats(1).avg_im = avg_im;
allstats(1).roi = roi;
allstats(1).thresholdFWHM = thresholdFWHM;
allstats(1).thresholdInv = thresholdInv;
allstats(1).CircVec = CircVec;
allstats(1).RadonMedFilt = RadonMedFilt;
allstats(1).vescenterx = windx;
allstats(1).vescentery = windy;
if current_cell == 1
    save([animal,'_PA',PA,'_',depth1,'_allstats.mat'],'allstats');
elseif current_cell == 2
    save([animal,'_PA',PA,'_',depth2,'_allstats.mat'],'allstats');
end

end

