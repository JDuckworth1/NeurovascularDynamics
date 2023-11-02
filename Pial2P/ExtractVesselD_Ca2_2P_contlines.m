%ExtractVesselD_Ca2__2P_contlines.m
% Load interleaved two photon images and (GCaMP channel and lumen Cy5.5
% channel). Calculate Ca df/f, integrated Cy5.5 intensity, and cross lines
% for FWH'M calculation. Perform background subtraction and average every 3
% frames of each channel.

%% 
clear; close all;
% data_folder = 'Y:\Data backup\20230613 JD221211F1';
data_folder = 'C:\Users\duckw\Desktop\FileTransfer'; %tiff stack location
savefolder = ''; %Folder to save processed data
animal = 'JD221223F2'; %For file naming

cd(data_folder);
files = dir('*001.tif'); %Get imaging tiff stack

file = 2;
namestr = extractBefore(files(file).name,'_00001.tif')
pix_um_x = 1; %1 um x resolution
pix_um_y = 4; %0.25 um y resolution
rate = input(['INPUT TRIAL FRAME RATE ',namestr,' (Hz) '])

%Load images
file_name = 'C:\Users\duckw\Desktop\FileTransfer\ROI1_00001.tif';
file_name_py = strrep(file_name,' ','**'); %Work around for spaces in file names
outvars = pyrunfile(sprintf("TiffDimDetector.py %s",file_name_py),'numims');
ims = str2double(outvars.char)

%Make vessel masks. Tried tiffreadVolume but was very slow compared to
%reading with python script
maskims = round(100*rate); %100 seconds * rate (images/s). We want multiple dilation cycles in mask
if mod(maskims,2) == 1 %If odd add 1 image for loading
    maskims = maskims + 1;
end

outvarsD = pyrunfile(sprintf("TiffImReader.py %s",file_name_py), 'imsout',r1 = int16(1), r2 = int16(maskims), r3 = int16(2));
Dmask = squeeze(mean(double(outvarsD),1)); %D mask is average of first approx 2000 images.

outvarsCa = pyrunfile(sprintf("TiffImReader.py %s",file_name_py), 'imsout',r1 = int16(0), r2 = int16(maskims - 1), r3 = int16(2));
Camask = squeeze(mean(double(outvarsCa),1));
% figure; subplot(2,1,1); imagesc(Dmask); daspect([1,4,1]); axis off; subplot(2,1,2); imagesc(Camask); daspect([1,4,1]); axis off;
fusemask = imfuse(Dmask,Camask); figure; imshow(fusemask); daspect([1,4,1]); 
cd(savefolder);
savefig('Ca_Lumen_imfust.fig');
%Save fuse image after resizing
fusemask_resize = imresize(fusemask,[size(fusemask,1),size(fusemask,2)*4],'nearest'); %Try not median filtering first
figure; imshow(fusemask_resize); daspect([1,1,1]); 
savefig('Ca_Lumen_imfust_resize.fig');
%Make combined mask 
Dmask = Dmask./max(Dmask(:)); Camask = Camask./max(Camask(:)); 
Combmask = Dmask + Camask;

if ~isfolder([animal,'_',namestr])
    mkdir([animal,'_',namestr]);
end
cd([animal,'_',namestr])

%Make mask from Ca + Lumen images
if ~isfile('Mask1.tif')
    figure; imshow(Combmask); daspect([1,4,1]);
    imwrite(Combmask,'Mask1.tif');
end
%Make mask from just lumen images.
if ~isfile('lumMask1.tif')
    figure; imshow(Dmask); daspect([1,4,1]);
    imwrite(Dmask,'lumMask1.tif');
end

%Plot std to help with mask creation if needed
stdmask = squeeze(std(double(outvarsD),[],1));
stdmask = stdmask./max(stdmask(:));
imwrite(stdmask,'stdmask2.tif');

%Use photoshop/GIMP to create binary mask from lumMask1.tif.
mask = logical(rgb2gray(imread('Mask.tif')));
camask = logical(rgb2gray(imread('CaMask.tif')));

%% Make vessel graph for Ca extraction
% mask1 is 0.25 x 1 um, mask is 0.25 x 0.25 resolution.

%For GCaMPintensity extraction make sure to use Ca-mask not lumen mask.
masksize = size(camask);
orig_masksize = size(Combmask);

opt.downsample_rate = 20; %Tunable parameter
opt.imdilate_disk_r = 1;
opt.min_kept_cc_num_pixel = 35;
opt.rm_neighbor_out_of_mask_Q = true;
opt.vis_mask_Q = true; % Show the cleaned up mask
[skel_neighbor_ind, ~] = fun_get_skeleton_neighbor_ind(camask, opt);
[skel_neighbor_ind_unique, tmp_unique_idx, ~] = unique(cat(2, skel_neighbor_ind{:}), 'stable');
recon_mask = false(size(camask));
recon_mask(skel_neighbor_ind_unique) = true;
num_skel_ind_unique = numel(skel_neighbor_ind_unique);
tmp_skel_label = repelem(1 : numel(skel_neighbor_ind), cellfun(@numel, skel_neighbor_ind));
skel_label_unique = tmp_skel_label(tmp_unique_idx);
toplot.skel_label = skel_label_unique;
toplot.mask_ind = skel_neighbor_ind_unique;

%Get images, average them and save matrix for wave extraction.
%rates range from 30.03 to 58.3 Hz. Average all data within 3 frames to get
%minimum rate of 10Hz.
%% JUST Load everything in, resize, and average appropriate number of images!!!
tic;
tmpimsCa = pyrunfile(sprintf("TiffImReader.py %s",file_name_py), 'imsout',r1 = int32(0), r2 = int32(round(ims)), r3 = int32(2));
Ca_data = double(tmpimsCa); clearvars tmpimsCa; %Ca_data = shiftdim(Ca_data,1); %Do this after background subtraction
toc
%% AVERAGE EVERY N frames to increase image quality while reducing frame rate

figure; imagesc(squeeze(Ca_data(100,:,:))); daspect([1,4,1]); title('Select rectangle for noise quantification');
[noisecol,noiserow] = ginput(2);
Ca_noise = Ca_data(:,round(noiserow(1)):round(noiserow(2)),round(noisecol(1)):round(noisecol(2)));
figure; histogram(Ca_noise); title('Choose Bounsd for gaussian fit');%Fit gaussian to noise
Ca_noise = Ca_noise(:);
[gausscol, ~] = ginput(2); 
Ca_noise_tofit = Ca_noise(and(Ca_noise > gausscol(1),Ca_noise < gausscol(2)));
pd = fitdist(Ca_noise_tofit,'Normal');
noise_mean = pd.mu %Noise level to be applied as offset manually.
clearvars Ca_noise_tofit

Ca_data = Ca_data - noise_mean; %Now noise is centered around zero (no photons = 0 intensity)
clearvars tmpimsD; Ca_data = shiftdim(Ca_data,1); %Make minimum value zero in Ca_data (int16)


FrameAvg = 3;
NewRate = rate/FrameAvg;
%First clip ending frames to make frames a multiple of FrameAvg
while mod(size(Ca_data,3),FrameAvg) ~= 0
    Ca_data(:,:,end) = [];
end
%Average Calcium data
Ca_data_tmp = reshape(Ca_data,[orig_masksize(1)*orig_masksize(2),size(Ca_data,3)])'; %collapse space dimension, now time x space
Ca_data_tmp = reshape(Ca_data_tmp,[FrameAvg,size(Ca_data,3)/FrameAvg*orig_masksize(1)*orig_masksize(2)]); %Each column is (FrameAvg) frames to be averaged
Ca_data_tmp = mean(Ca_data_tmp,1); %Take mean of every (FrameAvg) frames
Ca_data_tmp = reshape(Ca_data_tmp,[size(Ca_data,3)/FrameAvg,orig_masksize(1)*orig_masksize(2)])'; %Space (pixels) x time
Ca_data_tmp = reshape(Ca_data_tmp,[size(Ca_data,1),size(Ca_data,2),size(Ca_data,3)/FrameAvg]);
Ca_data = Ca_data_tmp; clearvars Ca_data_tmp;

%Convert isotropic skel_label to non-isotropic mask (0.25 x 1)um
map = zeros(size(camask));
map(toplot.mask_ind) = toplot.skel_label;
%Get set of big map indices for each small mask pixel. Set small mask pixel
%equal to mode of big map skel labels.
map_noni = zeros(size(Combmask));
for i = 1:size(Combmask,1)
    for j = 1:size(Combmask,2)
        mapcol = 4*(j-1)+1;
        map_noni(i,j) = mode(map(i,mapcol:mapcol + 3));
    end
end
%Make new skel_neighbor_ind for wave
% skel_neighbor_ind_noni = skel_neighbor_ind;
skel_neighbor_ind_noni = cell(1,1);
for i = 1:max(map_noni(:)) %Number of segments
    skel_neighbor_ind_noni{i,1} = find(map_noni == i)';
end

%CALCUATE Calcium dF/F from aveaged images
wave = fun_get_skeleton_neighbor_stat_from_image_stack(Ca_data, skel_neighbor_ind_noni, 'mean'); %Changed to mean for high mag imaging (most pixels have negligable change in fluorescence)
wave = bsxfun(@rdivide ,bsxfun(@minus, wave, mean(wave, 2)),mean(wave, 2));%calculate dF/F
time = 1:size(wave,2); time = time/NewRate; figure; plot(time,mean(wave,'omitnan')); xlabel('Time (s)','Interpreter','latex'); ylabel('GCaMP dF/F','Interpreter','latex');
%Save GCaMP data
h5create([animal,'_',namestr,'_wave.h5'],'/wave',[size(wave,1), size(wave,2)],'Datatype','double');
h5write([animal,'_',namestr,'_wave.h5'],'/wave', wave, [1 1],[size(wave,1), size(wave,2)]);
clearvars Ca_data

%% LOAD DIAMETER IMAGES AND AVERAGE THEM
tmpimsD = pyrunfile(sprintf("TiffImReader.py %s",file_name_py), 'imsout',r1 = int32(1), r2 = int32(round(ims+1)), r3 = int32(2));
D_data = double(tmpimsD);

figure; imagesc(squeeze(D_data(100,:,:))); daspect([1,4,1]); title('Select rectangle for noise quantification');
[noisecol,noiserow] = ginput(2);
D_noise = D_data(:,round(noiserow(1)):round(noiserow(2)),round(noisecol(1)):round(noisecol(2)));
figure; histogram(D_noise); title('Choose Bounsd for gaussian fit');%Fit gaussian to noise
D_noise = D_noise(:);
[gausscol, ~] = ginput(2); 
D_noise_tofit = D_noise(and(D_noise > gausscol(1),D_noise < gausscol(2)));
pd = fitdist(D_noise_tofit,'Normal');
noise_mean = pd.mu %Noise level to be applied as offset manually.
clearvars D_noise_tofit

D_data = D_data - noise_mean; %Now noise is centered around zero (no photons = 0 intensity)
clearvars tmpimsD; D_data = shiftdim(D_data,1); %Make minimum value zero in D_data (int16)
while mod(size(D_data,3),FrameAvg) ~= 0
    D_data(:,:,end) = [];
end
D_data_size = size(D_data);
% Average Lumen Data
D_data_tmp = reshape(D_data,[orig_masksize(1)*orig_masksize(2),D_data_size(3)])'; %collapse space dimension, now time x space
clearvars D_data
D_data_tmp = reshape(D_data_tmp,[FrameAvg,D_data_size(3)/FrameAvg*orig_masksize(1)*orig_masksize(2)]); %Each column is (FrameAvg) frames to be averaged
D_data_tmp = mean(D_data_tmp,1); %Take mean of every (FrameAvg) frames
D_data_tmp = reshape(D_data_tmp,[D_data_size(3)/FrameAvg,orig_masksize(1)*orig_masksize(2)])'; %Space (pixels) x time
D_data_tmp = reshape(D_data_tmp,[D_data_size(1),D_data_size(2),D_data_size(3)/FrameAvg]);
D_data = D_data_tmp; clearvars D_data_tmp;

%% USE CA SEGMENTS TO GET SUM OF INTENSITY OVER TIME

Dwave = fun_get_skeleton_neighbor_stat_from_image_stack(D_data, skel_neighbor_ind_noni, 'sum');
%Save Dwave
h5create([animal,'_',namestr,'_lumint.h5'],'/lumint',[size(Dwave,1), size(Dwave,2)],'Datatype','double');
h5write([animal,'_',namestr,'_lumint.h5'],'/lumint', Dwave, [1 1],[size(Dwave,1), size(Dwave,2)]);

%Calc cross correlation with gcamp signal
Dint_wave = mean(Dwave); Dint_wave = Dint_wave - mean(Dint_wave);
Caint_wave = mean(wave); Caint_wave = Caint_wave - mean(Caint_wave);
[corr,lags] = xcorr(Dint_wave,Caint_wave,'Normalized');
figure; plot(lags/NewRate,corr);

%% GET FWHM AT EACH SEGMENT, USE 10 CROSS LINES AVERAGE AT EACH SEGMENT
%Calculate cross lines in the equal pixel-size image, use
%imresize(,'nearest') to resize raw data and calculate profiles.
D_data = reshape(D_data,[size(D_data,1)*size(D_data,2),size(D_data,3)]);

h5create([animal,'_',namestr,'_diamdat.h5'],'/frame',[size(D_data,1), size(D_data,2)],'Datatype','double');
h5write([animal,'_',namestr,'_diamdat.h5'],'/frame', D_data, [1 1],[size(D_data,1), size(D_data,2)]);
%% Make cross lines for later diameter FWHM extraction.
% Use skeleton from vsl_graph, this is the common skeleton, everything is
% 1-1.
[~, ~, vsl_graph] = fun_get_skeleton_neighbor_ind(mask, opt);
[skel_neighbor_ind, ~] = fun_get_skeleton_neighbor_ind(mask, opt);
linkinds_tmp = cell2mat(vsl_graph.link.cc_ind);
maskskel = zeros(size(mask));
maskskel(linkinds_tmp) = 1;

%Now get radii
maskdist = bwdist(~mask);

[noderow,nodecol] = ind2sub(size(mask),vsl_graph.node.pos_ind);
figure; imagesc(maskskel); daspect([1,1,1]); f1 = gcf; f1.Position = [686 130 1021 957]; hold on;
plot(nodecol,noderow,'ro')
keepselecting = true;
counter = 1;
while keepselecting
    [x(counter),y(counter)] = ginput(1); 
    counter = counter + 1;
    keepselecting = input('Keep Selecting Links? 1 yes 0 no ');
end
%Find link that each point belongs to
islink = zeros(length(vsl_graph.link.cc_ind),1);
for i = 1:length(x)
    ptind = sub2ind(size(mask),round(y(i)),round(x(i)));
    for j = 1:length(vsl_graph.link.cc_ind)
        tmpinds = vsl_graph.link.cc_ind{j};
        %Dilate tmpinds to expand matching area
        tmpmap = zeros(size(mask)); tmpmap(tmpinds) = 1;
        tmpinds = find(imdilate(tmpmap,strel('disk',5)));
        if any(ismember(tmpinds,ptind))
            islink(j) = 1;
        end
    end
    i
end
linktodel = ~islink;
linkindcells = vsl_graph.link.cc_ind;
linkindcells(linktodel) = [];
linkinds = cell2mat(linkindcells); %Only inds of links to keep

keptmask = double(mask);
keptmask(linkinds) = 2;
figure; imagesc(keptmask); daspect([1,1,1]); title('Click first skeleton pixel')
[x_first,y_first] = ginput(1); 

%Downsample and calculate cross lines
%Instead of making continuous lines, put 10 lines at each segment (GCaMP
%and intensity data calculated at each segment already).

dmap = bwdistgeodesic(mask,round(x_first),round(y_first));
%Sort skeleton by increasing dmap values
linkinddistvals = dmap(linkinds);
[~,sortdist] = sort(linkinddistvals); 
sortlink = linkinds(sortdist); %link pixels put in order by distance along vessel

tmpmap = zeros(size(mask));
tmpmap(sortlink) = 1;

%Get link pixels to put crossline through. USE ONE SET OF LINES AT EACH SEG

numlines = length(skel_neighbor_ind); % #segments, Draw lines for each segment
crosslines = 10;
skelinds = linkinds; %Use the remaining links as our skeleton.
[skelrow,skelcol] = ind2sub(size(mask),skelinds);
%Get mean location & skel pixel for each segment
seg_skel = zeros(numlines,1);
for i = 1:numlines
    [current_row,current_col] = ind2sub(size(mask),skel_neighbor_ind{i});
    meanrow = mean(current_row);
    meancol = mean(current_col);

    meantoskel_dist = ((meanrow - skelrow).^2 + (meancol - skelcol).^2).^(1/2);
    [~,seg_skel(i)] = min(meantoskel_dist); %Minimum distance to each skeleton point. Gives skeleton for each segment.
end

%Delete seg_skel entries if they don't appear in sortlink
seg_skel_inds = skelinds(seg_skel);
[Lia,Locb] = ismember(seg_skel_inds,sortlink); %All seg_skel_inds should be member of sortlink.

seg_skel(~Lia) = []; %Now only inds in analyzed branch remain
seg_skel_inds(~Lia) = [];

%Put seg_skel values in order
Locb(~Lia) = [];
[~,Locb_sortinds] = sort(Locb);
seg_skel = seg_skel(Locb_sortinds); %Same order as sortlink.
seg_skel_inds = seg_skel_inds(Locb_sortinds);

%Get cross slopes for every ordered skeleton pixel. Then group by segment
xplot = 1:size(mask,2);
%Iterate over downsampled link pixels
[sortlink_crossrow,sortlink_crosscol] = ind2sub(size(mask),sortlink);
slopemean_num = 5;
midslope = nan(length(sortlink),1);
for i = 1:length(sortlink)
    currentind = sortlink(i); %Current link index to put crossline
    [current_row,current_col] = ind2sub(size(mask),currentind);
    if (i - slopemean_num) > 0 && (i + slopemean_num) < length(sortlink) %we can calculate from surrounding points
        meanrow_next = mean(sortlink_crossrow((i+1):(i+slopemean_num)));
        meanrow_prev = mean(sortlink_crossrow((i-slopemean_num):(i-1)));
        meancol_next = mean(sortlink_crosscol((i+1):(i+slopemean_num)));
        meancol_prev = mean(sortlink_crosscol((i-slopemean_num):(i-1)));
        midslope(i) = (meanrow_next-meanrow_prev)/(meancol_next-meancol_prev);
    end
end
crossslope = -1./midslope;
maskradii = maskdist(sortlink);

crosslengths_radii = 2.5

crosslengths = crosslengths_radii*maskradii;

counter = 1;
%Now we have slope for every skeleton pixel. Get ind cell array by
%collecting groups of (crosslines) at every segment pixel.
for i = 1:length(seg_skel_inds)
    skelpixtmp = seg_skel_inds(i);
    skelpix_loc = find(sortlink == skelpixtmp); %Where in the whole skeleton pixel list this segment's pixel falls
    %Get surrounding skeleton pixels
    if skelpix_loc < length(sortlink) - (crosslines-1) && skelpix_loc > (crosslines-1) %Make sure there is enough room for additional cross lines
        skelpix_surround = sortlink((skelpix_loc - (crosslines-1)):2:(skelpix_loc + (crosslines-1))); %Get surrounding skeleton pixels (skelpix_loc is the center)
        [skelpix_surround_row,skelpix_surround_col] = ind2sub(size(mask),skelpix_surround); %(x,y) coords of surrounding skeleton pixels

        singlevecs_row = skelpix_surround_row(2:end) - skelpix_surround_row(1:end-1);
        singlevecs_col = skelpix_surround_col(2:end) - skelpix_surround_col(1:end-1);
        meanvec = [mean(singlevecs_col),mean(singlevecs_row)]; %Col, row (x,y)
        meanvecslope = meanvec(2)/meanvec(1);
        meancrossslope = -1/meanvecslope;

        crossslopes = repmat(meancrossslope,[crosslines,1]);
        crosslength = crosslengths((skelpix_loc - (crosslines-1)):2:(skelpix_loc + (crosslines-1)));
        meancrosslength = mean(crosslength);

        for j = 1:length(crosslength)
            if abs(crossslopes(j)) ~= Inf
                ycross_eq = crossslopes(j).*(xplot - skelpix_surround_col(j)) + skelpix_surround_row(j);
                [xcross_subs,ycross_subs] = bresenham(xplot(1),ycross_eq(1),xplot(end),ycross_eq(end));
            else %Infinite slope -> vertical line
                ycross_subs = 1:size(mask,1);
                xcross_subs = repmat(skelpix_surround_col(j),[1,length(ycross_subs)]);
            end
            crossdists = ((ycross_subs - skelpix_surround_row(j)).^2 + (xcross_subs - skelpix_surround_col(j)).^2).^(1/2);
            keepcross = crossdists < meancrosslength;
            %Check if any pixels are outside of image
            xcross_subs = xcross_subs(keepcross); ycross_subs = ycross_subs(keepcross);
            if any(xcross_subs > size(mask,2)) || any(ycross_subs > size(mask,1)) || any(xcross_subs < 1) || any(ycross_subs < 1)
                linetodel(counter) = 1;
                counter = counter + 1;
            else
                linetodel(counter) = 0;
                inds{counter} = sub2ind(size(mask),ycross_subs,xcross_subs);
                counter = counter + 1;
            end
        end
    end
end
inds = inds(~cellfun('isempty',inds));
howuneven = mod(length(inds),crosslines);
if howuneven ~= 0
    cellsz = cell2mat(cellfun(@length,inds,'uni',false));
    if cellsz(end-(howuneven)) ~= cellsz(end-(howuneven-1)) && howuneven ~= 0
        inds((end-(howuneven - 1)):end) = [];
    elseif cellsz(howuneven+1) ~= cellsz(howuneven) && howuneven ~= 0
        inds(1:howuneven) = [];
    end
end

%Keep only unique inds
%Make array with each row is one set of lines. Then use 
indscheck = zeros(max(cell2mat(cellfun(@length,inds,'uni',false))),length(inds)/crosslines);
group = 1;
for i = 1:length(inds)
    if mod(i,crosslines) == 1 %Just look at the first line in each group
        indscheck(1:length(inds{1,i}),group) = inds{1,i};
        group = group + 1;
    end
end
indscheck = indscheck';
[~,indsia,~] = unique(indscheck,'rows'); %Inds
indstokeep = (indsia * crosslines) - (0:(crosslines-1));
indstokeep = sort(indstokeep(:));
inds = inds(indstokeep); %Only keep unique line groups

maskplot = double(logical(rgb2gray(imread('Mask.tif'))));
for i=1:length(inds)
    maskplot(inds{i}) = 2;
end
maskplot(skelinds) = 2;
figure; imagesc(maskplot); daspect([1,1,1]);

save('seg_inds.mat','inds');

% SAVE PLOTTING INFO AND EXPERIMENT INFO
vsl_graph.link.linktodel = linktodel;

toplot.sortlink = sortlink;
toplot.seg_skel_inds = seg_skel_inds;
toplot.vsl_graph = vsl_graph;
toplot.rate = rate;
toplot.NewRate = NewRate;
toplot.crosslengths_radii = crosslengths_radii;
toplot.crosslines = crosslines;
toplot.frameavg = FrameAvg;
toplot.manualoffset = noise_mean;
save('toplot.mat','toplot')


