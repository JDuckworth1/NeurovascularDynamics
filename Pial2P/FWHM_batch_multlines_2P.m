% FWHM_batch_multlines_2P.m

% Use previously calculated cross lines to calculate the FW(scalefactor)Max(t) at every
% location along a given pial vessel. Written for use with parpool.

%%
roiname = 'ROI1'
animal = 'JD221223F2'
filepath = [animal,'/',roiname];
cd(filepath)

tmpdiam = struct(); %For saving temporary extracted diameters
load('seg_inds.mat');
origmask = imread('Mask_origsize.tif');
mask = imread('Mask.tif');
cd(filepath);
d_data = h5read([animal,'_',roiname,'_diamdat.h5'],'/frame');
cd(filepath)
d_data = reshape(d_data,[size(origmask,1),size(origmask,2),size(d_data,2)]); %X x Y x Time
disp('Done Loading')

meanimage = mean(d_data,3);

%%
crosslines = 10;
numlinegroups = length(inds)/crosslines;

linesiter = 1:crosslines:numlinegroups*crosslines;
diam = zeros(size(d_data,3),numlinegroups);
pix_mm = 4000;

cellsizes = cellfun(@length,inds,'UniformOutput',true);
maxcellsize = max(cellsizes);
indmat = NaN(maxcellsize,length(inds));
for i = 1:length(inds)
    indmat(1:length(inds{i}),i) = inds{i};
end
indmat_todel = isnan(indmat);
indmat(indmat_todel) = 1; %Temporarily, get the first pixel of the image, then delete later.

pix_mm_lines = zeros(numlinegroups,1);
for i = 1:numlinegroups
    linetmp = indmat(:,linesiter(i));
    linetmp = linetmp(~indmat_todel(:,linesiter(i)));
    [tmprow,tmpcol] = ind2sub(size(mask(:,:,1)),linetmp);
    pix_dist_line = sqrt((tmprow(end)-tmprow(1))^2 + (tmpcol(end)-tmpcol(1))^2); %Units: pixels
    pix_mm_lines(i) = length(linetmp)/((1/pix_mm)*pix_dist_line); %pixels in line per mm
end
pix_mm_lines = pix_mm_lines';

%Analysis parameters
sf = 0.3;
MedFiltVal = 9; %For image
ProfMedFiltVal = 6; %For average profile

pools = 48;
parpool(pools);

tic
parfor j = 1:size(d_data,3) %Iterate over frames

    numlines = numlinegroups;
    imtmp = imresize(d_data(:,:,j),[size(mask,1),size(mask,2)],'nearest');
    imtmp = medfilt2(imtmp,[MedFiltVal,MedFiltVal]);
    imtmp = rescale(imtmp);

    counter = 1;
    proftmp = nan();
    mincellsize = nan(numlinegroups,1);
    for line = linesiter
        cellsizes_tmp = cellsizes(line:line+(crosslines-1));
        mincellsize(counter) = min(cellsizes_tmp);
        proftmp(1:mincellsize(counter),counter) = mean(double(imtmp(indmat(1:mincellsize(counter),line:line+(crosslines-1)))),2); %Average of 10 line profiles
        %             proftmp(1:mincellsize(counter),counter) = median(double(imtmp(indmat(1:mincellsize(counter),line:line+(crosslines-1)))),2); %Try this for trials with noisy estimates 
        counter = counter + 1;
    end
    tonan = proftmp == 0;
    proftmp = medfilt1(proftmp,ProfMedFiltVal,'omitnan');
    proftmp(tonan) = nan;

    proftmp_flip = NaN(size(proftmp));
    for i = 1:numlines
        tmpcol = proftmp(:,i);
        proftmp_flip(1:mincellsize(i),i) = flip(tmpcol(~isnan(tmpcol)));
    end
    proftmp_flip(11:end,:) = []; %Use last 10 entries to estimate minimum
    profstart = proftmp(1:10,:);
    
    Minvaltest = [std(profstart,[],1);std(proftmp_flip,[],1)];
    [~,Minval_isflatprof] = min(Minvaltest,[],1); %1 if start of profile is less variable, 2 if end of profile is less variable.
    Minval_isprofstart = Minval_isflatprof == 1; %Logical, which lines to use start of profile for.
    Minval = zeros(1,size(proftmp_flip,2));
    Minval(Minval_isprofstart) = median(profstart(:,Minval_isprofstart));
    Minval(~logical(Minval_isprofstart)) = median(proftmp_flip(:,~logical(Minval_isprofstart)));

    Maxval = max(proftmp,[],1,'omitnan'); %Max pixel intensity of each profile
    HMval = (Maxval - Minval)*sf + Minval;

    intpfactor = 0.01;
    %Interpolate so that we have 100x more points than starting
    %(each profile stars with silghtly different # of points, here mm/pixel will be constant)
    intpsizes = (mincellsize - 1).*(1/intpfactor)+1; %Size of interpolated profiles
    sprofintp = nan(max(intpsizes),numlines);
    for i = 1:numlines
        sprofintp(1:intpsizes(i),i) = interp1(1:mincellsize(i),proftmp(1:mincellsize(i),i),1:intpfactor:mincellsize(i));
    end

    intppix_mm = pix_mm_lines/intpfactor; %Each interpolated pixel length is (intpfactor)*pixel length

    tol = 0.005; %This does not have to be large after profile interpolation. Need to make sure vessel midpoint dip isn't within threshold.
    HMval_repmat = repmat(HMval,[size(sprofintp,1),1]);
    prof_diff = abs(sprofintp - HMval_repmat);
    v = zeros(1,numlines);
    for i = 1:numlines
        small_diffs = find(prof_diff(:,i) < tol); %indices of profile distance to half-max less than tolerance

        [~,groupsep] = max(diff(small_diffs)); %last ind of group 1

        [~,group1ID] = min(abs(sprofintp(small_diffs(1:groupsep)) - HMval(i))); %Ind of closest group 1 point, within group 1
        group1ID = group1ID + small_diffs(1) - 1; %Ind of closest group 1 point, within total profile
        [~,group2ID] = min(abs(sprofintp(small_diffs(groupsep+1:end)) - HMval(i))); %Ind of closest group 2 point, within group 2
        group2ID = group2ID + small_diffs(groupsep + 1) - 1; %Ind of closest group 2 point, within total profile

        v(i) = (group2ID - group1ID)/intppix_mm(i);
    end
    diam(j,:) = v;
end
toc

%Save diameter and parameters used. 
tmpdiam.diam = diam; %Here beforepuff is col 1 puff is 2 afterpuff is 3
tmpdiam.sf = sf; %FW(scalefactor)Max 
tmpdiam.MedFiltVal = MedFiltVal;
tmpdiam.ProfMedFiltVal = ProfMedFiltVal;

disp('FWHM Done');
save('tmpdiam_2P.mat','tmpdiam')

