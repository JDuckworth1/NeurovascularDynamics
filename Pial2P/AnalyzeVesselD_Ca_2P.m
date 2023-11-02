%AnalyzeVesselD_Ca_2P.m
clear; clc; close all;

roiname = 'ROI1'
animal = 'JD221223F2'
data_folder = [animal,'\',roiname];
cd(data_folder);

%Load diameter FW(Scalefactor)M
diamfiles = dir('tmpdiam_2P.mat');
load(diamfiles(1).name); 
diamFWHM = tmpdiam.diam;
load('seg_inds.mat');
%Load calcium dF/F
wavefiles = dir('*wave.h5');
wave = h5read(wavefiles(1).name,'/wave');
%Load lumen sum of intensity data
lumenfiles = dir('*lumint.h5');
lumwave = h5read(lumenfiles(1).name,'/lumint');
%Load plotting and experiment details
toplotfiles = dir('*toplot.mat');
load(toplotfiles(1).name);
%Load mask
mask = logical(rgb2gray(imread('Mask.tif')));

clearvars -except data_folder diamFWHM wave lumwave toplot inds mask roiname
time = 1:size(lumwave,2);
time = time/toplot.NewRate;
time = (time - time(1)) + (1/(2*toplot.NewRate));

figure
plot(time,mean(diamFWHM,2)/max(mean(diamFWHM)),'k');
hold on
plot(time,mean(lumwave)/max(mean(lumwave)),'b');
legend({'Diam FWHM','Sum Lumen Intensity'})
xlabel('Time (s)','Interpreter','latex');
ylabel('Normalized units','Interpreter','latex');

figure
plot(time,mean(wave,'omitnan')/max(mean(wave,'omitnan')),'g');
hold on
plot(time,mean(lumwave)/max(mean(lumwave)),'b');
legend({'Ca dF/F','Sum Lumen Intensity'})
xlabel('Time (s)','Interpreter','latex');
ylabel('Normalized units','Interpreter','latex');


mask_inds = toplot.mask_ind;
found_skel_labels = zeros(length(inds)/toplot.crosslines,1);
for i = 1:length(inds)/toplot.crosslines %Iterate over lines
    currentlinenum = (i*toplot.crosslines) - (toplot.crosslines - 1);
    for j = 1:toplot.crosslines
        if j == 1
            singlineinds = inds{currentlinenum};
            if size(singlineinds,1) == 1
                singlineinds = singlineinds';
            end
            tmpinds = singlineinds;
        else
            singlineinds = inds{currentlinenum+j-1};
            if size(singlineinds,1) == 1
                singlineinds = singlineinds';
            end
            tmpinds = [tmpinds; singlineinds];
        end
    end
    [tmprow,tmpcol] = ind2sub(size(mask),tmpinds);
    meanrow = mean(tmprow);
    meancol = mean(tmpcol);
    meanind = sub2ind(size(mask),round(meanrow),round(meancol));

    %Check which skeleton seg this is in
    [~,locmaskind] = ismember(mask_inds,meanind);
    locmaskind = find(locmaskind);
    found_skel_labels(i) = toplot.skel_label(locmaskind); %Which wave row this line corresponds to.
end

[~,uniqueb,~] = unique(found_skel_labels);
repeatskel_labels = 1:size(diamFWHM,2);
repeatskel_labels(uniqueb) = nan;
repeatskel_labels(isnan(repeatskel_labels)) = [];
indstodel = zeros(length(inds),1);
%Replace repeats with one
for i = 1:length(repeatskel_labels)
    tmpinds = [repeatskel_labels(i)-1,repeatskel_labels(i)];
    FWHMtmpmat = [diamFWHM(:,repeatskel_labels(i)-1),diamFWHM(:,repeatskel_labels(i))];
    tmpstd = std(FWHMtmpmat,[],1);
    [~,stddelete] = max(tmpstd);
    %Delete diamFWHM and inds
    diamFWHM(:,tmpinds(stddelete)) = nan;
    found_skel_labels(tmpinds(stddelete)) = nan;

    currentlinenum = tmpinds(stddelete)*toplot.crosslines-(toplot.crosslines-1);
    indstodel(currentlinenum:(currentlinenum + toplot.crosslines - 1)) = 1;
end
isnanvec = isnan(diamFWHM(1,:));
diamFWHM(:,isnanvec) = [];
found_skel_labels(isnanvec) = [];
%delete repeated inds
inds(logical(indstodel)) = [];

%order found_skel_labels in increasing order
[found_skel_labs_sort,FSLInds] = sort(found_skel_labels);
%Sort diamFWHM
diamFWHM = diamFWHM(:,FSLInds);
%Don't need to sort lumwave or wave. They are ordered by increasing
%skeleton segment.

%Sort inds
linessort2 = FSLInds*toplot.crosslines;
linessort1 = linessort2 - 9;
linessort = [];
for i = 1:length(linessort1)
    linessort = [linessort,linessort1(i):linessort2(i)];
end
inds = inds(linessort);
%Delete rows from wave and lumwave
wavetokeep = zeros(size(wave,1),1);
wavetokeep(found_skel_labs_sort) = 1;
wave(~logical(wavetokeep),:) = [];
lumwave(~logical(wavetokeep),:) = [];

diamFWHM = diamFWHM'; %Now diamFWHM, wave, lumwave are (segments x samples), inds are ordered by skeleton segment as well.
if ~isequal(size(diamFWHM,1),size(wave,1),size(lumwave,1))
    disp('Data sizes arent equal'); pause();
end
clearvars -except data_folder diamFWHM wave lumwave toplot inds mask time found_skel_labs_sort wavetokeep indstodel roiname isnanvec


% Time x space analysis
%DOES LUMEN INTENSITY (AREA) SCALE AS FWHM (DIAMETER)^2?
diam = diamFWHM * 1000; %Change units to um
% EXCLUDE SEGMENTS WITH POOR DIAMETER EXTRACTION use previous function
% Can also exclude movement periods
% addpath('X:\DataAnalysis\VesCorrPhase\TxRed_GCaMP');
[todel,~,~,~,~,~,~] = preprocessdiam_VisCrosslines(diam',toplot.NewRate,inds,mask);

diam(todel,:) = [];
wave(todel,:) = [];
lumwave(todel,:) = [];
found_skel_labs_sort(todel) = []; %Tells us the skel_label for each wave row

% subtract mean from all data
wavemean = mean(wave,2);
lumwavemean = mean(lumwave,2);
diammean = mean(diam,2);
wavemeanmat = repmat(wavemean,[1,size(wave,2)]);
lumwavemeanmat = repmat(lumwavemean,[1,size(lumwave,2)]);
diammeanmat = repmat(diammean,[1,size(diam,2)]);
wave_zmean = wave - wavemeanmat;
lumwave_zmean = lumwave - lumwavemeanmat;
diam_zmean = diam - diammeanmat;
% clearvars diammeanmat wavemeanmat lumwavemeanmat

%Low pass filter all data
lowp = designfilt('lowpassiir','PassbandFrequency',0.3,...
    'StopbandFrequency',0.35,'StopbandAttenuation',10,'SampleRate',toplot.NewRate);
wave_lp = transpose(filtfilt(lowp,wave_zmean'));
lumwave_lp = transpose(filtfilt(lowp,lumwave_zmean'));
diam_lp = transpose(filtfilt(lowp,diam_zmean'));

wave_lp = wave_lp + wavemeanmat;
lumwave_lp = lumwave_lp + lumwavemeanmat;
diam_lp = diam_lp + diammeanmat;

%Calculate intensity - diameter relationship at each location
plotlumD = 0;
for i = 1:size(diam,1)
    xfit = log10(diam_lp(i,:));
    yfit = log10(lumwave_lp(i,:));
    meany = mean(yfit);
%     xfit = log10(diam(i,:));
%     yfit = log10(lumwave(i,:));
    fittmp = polyfit(xfit,yfit,1);
    I_DSlope(i) = fittmp(1);
    I_DConst(i) = fittmp(2);
        y_est = polyval(fittmp,xfit); hold on; plot(xfit,y_est,'r--');
    SSR = sum((yfit - y_est).^2);
    SST = sum((yfit - meany).^2);
    I_DR2(i) = 1-SSR/SST;
    if plotlumD == 1
        figure; scatter(xfit,yfit,'filled','MarkerFaceAlpha',0.05);
        title([num2str(i),' ',num2str(I_DSlope(i))]);
        ylim([min(min(log10(lumwave_lp(:,:)))),max(max(log10(lumwave_lp(:,:))))])
        xlim([min(min(log10(diam_lp(:,:)))),max(max(log10(diam_lp(:,:))))])
        %     ylim([min(min(log10(lumwave(:,:)))),max(max(log10(lumwave(:,:))))])
        %     xlim([min(min(log10(diam(:,:)))),max(max(log10(diam(:,:))))])
    end
end

% Calculate phase progressions

%Get distances from segment locations
% cellsizes = cellfun(@length,inds,'UniformOutput',true);
rate = toplot.NewRate;
% rate = (39.64/3);
addpath('X:\DataAnalysis\VesCorrPhase\TxRed_GCaMP')

%Iterate over segments
skelx = zeros(size(wave,1),1); skely = zeros(size(wave,1),1);
for i=1:size(wave,1)
    skelindsloc = toplot.skel_label == found_skel_labs_sort(i);
    skelmaskinds = toplot.mask_ind(skelindsloc);
    [skelrows,skelcols] = ind2sub(size(mask),skelmaskinds);

    skelx(i) = mean(skelcols(:),'omitnan');
    skely(i) = mean(skelrows(:),'omitnan');
end

%Calc Distances
pix_mm = 4000; %0.25um/pixel
dmap = bwdistgeodesic(mask,round(skelx(1)),round(skely(1)));
numpts = size(wave,1);
geodist = zeros(numpts,1);
for i=1:numpts
    geodist(i) = dmap(round(skely(i)),round(skelx(i)));
end
geodist = geodist/pix_mm; %Distances in mm


%Find vasomotion peaks before calculating phases
params.Fs = rate; %Interpolated rate is twice actual single depth rate
params.pad = 2;
params.fpass = [0 params.Fs/2]; %Hz, default is [0 Fs/2]
params.err   = [2 .05];
params.trialave = 0;
T = time(end);
BW = 0.02; %500s trial -> TBW =~ 10
if T < 300
    BW = 0.04; %Keep TBW ~= 10 for 250s trial
end
params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers
addpath(genpath('\chronux_2_12'))
[S,f] = mtspectrumc(diam_zmean',params);
[CaS,Caf] = mtspectrumc(wave_zmean',params);
%get Calcium vasomotion peak
figure
plot(Caf,log10(mean(CaS,2))); xlim([0 0.5]);
[wind,~] = ginput(2);
findf1 = round(wind(1),3);
findf2 = round(wind(2),3);
f1 = find(round(Caf,3)==round(findf1,3),1,'first'); %Hz lower bound
f2 = find(round(Caf,3)==round(findf2,3),1,'first'); %Hz upper bound
rmpath(genpath('\chronux_2_12'))
[maxes,f_peaks] = findpeaks(CaS(f1:f2),Caf(f1:f2));
maxloc = find(maxes==max(maxes));
f_peak_Ca = f_peaks(maxloc)
freq_findx_Ca = max(find(round(Caf,3)==round(f_peak_Ca,3)));
addpath(genpath('\chronux_2_12'))
Ca_fv_findfv = [f_peak_Ca,freq_findx_Ca];

%Same for diameter vasomotion peak
figure
plot(f,log10(mean(S,2))); xlim([0 0.5]);
[wind,~] = ginput(2);
findf1 = round(wind(1),3);
findf2 = round(wind(2),3);
f1 = find(round(f,3)==round(findf1,3),1,'first'); %Hz lower bound
f2 = find(round(f,3)==round(findf2,3),1,'first'); %Hz upper bound
rmpath(genpath('\chronux_2_12'))
[maxes,f_peaks] = findpeaks(S(f1:f2),f(f1:f2));
maxloc = find(maxes==max(maxes));
f_peak_D = f_peaks(maxloc)
freq_findx_D = max(find(round(f,3)==round(f_peak_D,3)));
addpath(genpath('\chronux_2_12'))
D_fv_findfv = [f_peak_D,freq_findx_D];


%Diameter FWHM
% Calculate phase progressions
[tapers1,pad1,Fs,fpass1,err1,trialave1]=getparams(params);
N1 = size(wave,2);
nfft1=max(2^(nextpow2(N1)+pad1),N1);
[f1,findx1]=getfgrid(Fs,nfft1,fpass1);
tapers1=dpsschk(tapers1,N1,Fs);
wave_forcalc = diam_zmean';
for i=1:numpts %Iterate over segments
    if i==1 %If we're on the first segment
        wave0 = wave_forcalc(:,i);
        Jneu0 = mtfftc(wave0,tapers1,nfft1,Fs);
        Jneuplot0 = Jneu0(findx1,:,:);
        clear Jves0

        wave1 = wave_forcalc(:,i);
        Jneu1 = mtfftc(wave1,tapers1,nfft1,Fs);
        Jneuplot1 = Jneu1(findx1,:,:);
        clear Jves1

        S12_Ca = squeeze(mean(conj(Jneuplot1(freq_findx_Ca,:)).*Jneuplot0(freq_findx_Ca,:),2));
        phi_mat_D(i) = angle(S12_Ca);
        S12_Ca_fD = squeeze(mean(conj(Jneuplot1(freq_findx_D,:)).*Jneuplot0(freq_findx_D,:),2));
        phi_mat_D_fD(i) = angle(S12_Ca_fD);
    else
        wave1 = wave_forcalc(:,i);
        Jneu1 = mtfftc(wave1,tapers1,nfft1,Fs);
        Jneuplot1 = Jneu1(findx1,:,:);
        clear Jves1

        S12_Ca = squeeze(mean(conj(Jneuplot1(freq_findx_Ca,:)).*Jneuplot0(freq_findx_Ca,:),2));
        phi_mat_D(i) = angle(S12_Ca);
        S12_Ca_fD = squeeze(mean(conj(Jneuplot1(freq_findx_D,:)).*Jneuplot0(freq_findx_D,:),2));
        phi_mat_D_fD(i) = angle(S12_Ca_fD);
    end
end
% figure
% scatter(geodist,phi_mat_D,'filled');
% xlabel('Distance (mm)','Interpreter','latex');
% ylabel('Phase (rad','Interpreter','latex');
% title('Diameter FWHM phase progression','Interpreter','latex');

%Sum lumen intensity
wave_forcalc = lumwave';
for i=1:numpts %Iterate over segments
    if i==1 %If we're on the first segment
        wave0 = wave_forcalc(:,i);
        Jneu0 = mtfftc(wave0,tapers1,nfft1,Fs);
        Jneuplot0 = Jneu0(findx1,:,:);
        clear Jves0

        wave1 = wave_forcalc(:,i);
        Jneu1 = mtfftc(wave1,tapers1,nfft1,Fs);
        Jneuplot1 = Jneu1(findx1,:,:);
        clear Jves1

        S12_Ca = squeeze(mean(conj(Jneuplot1(freq_findx_Ca,:)).*Jneuplot0(freq_findx_Ca,:),2));
        phi_mat_lumD(i) = angle(S12_Ca);
        S12_Ca_fD = squeeze(mean(conj(Jneuplot1(freq_findx_D,:)).*Jneuplot0(freq_findx_D,:),2));
        phi_mat_lumD_fD(i) = angle(S12_Ca_fD);
    else
        wave1 = wave_forcalc(:,i);
        Jneu1 = mtfftc(wave1,tapers1,nfft1,Fs);
        Jneuplot1 = Jneu1(findx1,:,:);
        clear Jves1

        S12_Ca = squeeze(mean(conj(Jneuplot1(freq_findx_Ca,:)).*Jneuplot0(freq_findx_Ca,:),2));
        phi_mat_lumD(i) = angle(S12_Ca);
        S12_Ca_fD = squeeze(mean(conj(Jneuplot1(freq_findx_D,:)).*Jneuplot0(freq_findx_D,:),2));
        phi_mat_lumD_fD(i) = angle(S12_Ca_fD);
    end
end
% figure
% scatter(geodist,phi_mat_lumD,'filled');
% xlabel('Distance (mm)','Interpreter','latex');
% ylabel('Phase (rad','Interpreter','latex');
% title('Diameter Integrated Intensity phase progression','Interpreter','latex');

%Calcium
wave_forcalc = wave';
for i=1:numpts %Iterate over segments %Want reference to be non-conj so that if it leads then phase is positive.
    if i==1 %If we're on the first segment
        wave0 = wave_forcalc(:,i);
        Jneu0 = mtfftc(wave0,tapers1,nfft1,Fs);
        Jneuplot0 = Jneu0(findx1,:,:);
        clear Jves0

        wave1 = wave_forcalc(:,i);
        Jneu1 = mtfftc(wave1,tapers1,nfft1,Fs);
        Jneuplot1 = Jneu1(findx1,:,:);
        clear Jves1

        S12_Ca = squeeze(mean(conj(Jneuplot1(freq_findx_Ca,:)).*Jneuplot0(freq_findx_Ca,:),2));
        phi_mat_Ca(i) = angle(S12_Ca);
        S12_Ca_fD = squeeze(mean(conj(Jneuplot1(freq_findx_D,:)).*Jneuplot0(freq_findx_D,:),2));
        phi_mat_Ca_fD(i) = angle(S12_Ca_fD);
    else
        wave1 = wave_forcalc(:,i);
        Jneu1 = mtfftc(wave1,tapers1,nfft1,Fs);
        Jneuplot1 = Jneu1(findx1,:,:);
        clear Jves1

        S12_Ca = squeeze(mean(conj(Jneuplot1(freq_findx_Ca,:)).*Jneuplot0(freq_findx_Ca,:),2));
        phi_mat_Ca(i) = angle(S12_Ca);
        S12_Ca_fD = squeeze(mean(conj(Jneuplot1(freq_findx_D,:)).*Jneuplot0(freq_findx_D,:),2));
        phi_mat_Ca_fD(i) = angle(S12_Ca_fD);
    end

    tmpca = wave_lp(1,:) - wavemean(1);
    tmpca1 = wave_lp(i,:) - wavemean(i);
    [corrtmp,lagstmp] = xcorr(tmpca1,tmpca,'normalized'); %Max lag positive when intensity lags ca2+
    [maxcorr,maxlag(i)] = max(corrtmp);
end
figure
scatter(geodist,lagstmp(maxlag)/rate,'filled');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Calcium lag progression','Interpreter','latex');
% ylabel('Phase (rad','Interpreter','latex');
% title('Diameter Calcium cphase progression','Interpreter','latex');

%Calculate slope from fit and errors, R2
[phasevec_D] = MakePhasevec_Cy55GCaMP(geodist,phi_mat_D',f_peak_D);
[phasevec_D_fD] = MakePhasevec_Cy55GCaMP(geodist,phi_mat_D_fD',f_peak_Ca);
[phasevec_I] = MakePhasevec_Cy55GCaMP(geodist,phi_mat_lumD',f_peak_D);
[phasevec_I_fD] = MakePhasevec_Cy55GCaMP(geodist,phi_mat_lumD_fD',f_peak_Ca);
[phasevec_Ca] = MakePhasevec_Cy55GCaMP(geodist,phi_mat_Ca',f_peak_D);
[phasevec_Ca_fD] = MakePhasevec_Cy55GCaMP(geodist,phi_mat_Ca_fD',f_peak_Ca);

% Calculate relation between calcium and constriction at each location

%Integrated intensity and Calcium
for i = 1:numpts
    tmpintensity = lumwave_lp(i,:) - lumwavemean(i);
    tmpca = wave_lp(i,:) - wavemean(i);
    [corrtmp,lagstmp] = xcorr(tmpintensity,tmpca,'normalized'); %Max lag positive when intensity lags ca2+

    [maxcorr,maxlag] = min(corrtmp);
    I_Ca_lagvec(i) = lagstmp(maxlag)/rate;
    I_Ca_corrvec(i) = maxcorr;
end

for i = 1:numpts
    tmpD = diam_lp(i,:) - diammean(i);
    tmpca = wave_lp(i,:) - wavemean(i);
    [corrtmp,lagstmp] = xcorr(tmpD,tmpca,'normalized'); %Max lag positive when intensity lags ca2+

    [maxcorr,maxlag] = min(corrtmp);
    D_Ca_lagvec(i) = lagstmp(maxlag)/rate;
    D_Ca_corrvec(i) = maxcorr;
end


% Save results
resultsmat(1).diam = diam;
resultsmat(1).diammean = diammean;
resultsmat(1).lumwave = lumwave;
resultsmat(1).lumwavemean = lumwavemean;
resultsmat(1).wave = wave;
resultsmat(1).wavemean = wavemean;
resultsmat(1).rate = rate;
resultsmat(1).time = time;
resultsmat(1).geodist = geodist;
resultsmat(1).phi_mat_D = phi_mat_D;
resultsmat(1).phi_mat_lumD = phi_mat_lumD;
resultsmat(1).phi_mat_Ca = phi_mat_Ca;
resultsmat(1).phi_mat_D_fD = phi_mat_D_fD;
resultsmat(1).phi_mat_lumD_fD = phi_mat_lumD_fD;
resultsmat(1).phi_mat_Ca_fD = phi_mat_Ca_fD;
resultsmat(1).I_Ca_lagvec = I_Ca_lagvec;
resultsmat(1).I_Ca_corrvec = I_Ca_corrvec;
resultsmat(1).D_Ca_lagvec = D_Ca_lagvec;
resultsmat(1).D_Ca_corrvec = D_Ca_corrvec;
resultsmat(1).phasevec_D = phasevec_D;
resultsmat(1).phasevec_I = phasevec_I;
resultsmat(1).phasevec_Ca = phasevec_Ca;
resultsmat(1).phasevec_D_fD = phasevec_D_fD;
resultsmat(1).phasevec_I_fD = phasevec_I_fD;
resultsmat(1).phasevec_Ca_fD = phasevec_Ca_fD;
resultsmat(1).lpfilt = lowp;
resultsmat(1).I_DSlope = I_DSlope;
resultsmat(1).I_DConst = I_DConst;
resultsmat(1).I_DR2 = I_DR2;
resultsmat(1).todel = todel;
resultsmat(1).indstodel = indstodel;
resultsmat(1).isnanvec = isnanvec;

save('resultsmat.mat','resultsmat');

%% Plot results
colormat = zeros(numpts,3); colormat(:,1) = 0; colormat(:,3) = 1; colormat(:,2) = flip((1/numpts):(1/numpts):1);
coloriter = 1;
%Plot locations on vessel image
h1 = openfig('Ca_Lumen_imfust.fig','reuse'); 
ax1 = gca;
for i = 1:length(inds)/10
    if todel(i) == 0
        hold on;
        linenum = i*10-5;
        ind1 = inds{1,linenum}(1);
        ind2 = inds{1,linenum}(end);
        [row1,col1] = ind2sub(size(mask),ind1);
        [row2,col2] = ind2sub(size(mask),ind2);
        line([col1,col2]/4,[row1,row2],'Color',colormat(coloriter,:))
        coloriter = coloriter + 1;
    end
end
fig1 = get(ax1,'children');

fig=figure('units','inches','outerposition',[0 0 8.5 11]); hold on;
s1 = subplot(5,2,[1 2]);
copyobj(fig1,s1);
ylim([0 size(mask,1)]); xlim([0 size(mask,2)/4]);
daspect([1,4,1]);
axis off

subplot(5,2,[3 4])
yyaxis('left')
plot(time,mean(diam,1),'k');
ylabel('Diam (um)','Interpreter','latex','Color','k');
set(gca,'ycolor','k') 
yyaxis right
plot(time,mean(wave),'g');
ylabel('Ca2+ dF/F','Interpreter','latex','Color','g');
set(gca,'ycolor',[0.4660 0.6740 0.1880]) 
set(gca, 'YDir','reverse')
legend({'Diam FWHM','Calcium dF/F'})
xlabel('Time (s)','Interpreter','latex');

%Plot phase progressions
subplot(5,2,[5 6]);
scatter(geodist,phi_mat_D,[],colormat,'filled','k');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Phase (rad)','Interpreter','latex');
% title(['Diameter FWHM phase progression, v = ',num2str(round(resultsmat(1).phasevec_D(7),2)),' mm/s, $R^2$ = ',num2str(round(resultsmat(1).phasevec_D(2)^2,2))],'Interpreter','latex');
xplot = 0:0.001:max(geodist); yplot = resultsmat(1).phasevec_D(1)*xplot + resultsmat(1).phasevec_D(8);
hold on; % plot(xplot,yplot,'k');
%Intensity sum
% subplot(5,2,[5 6]);
scatter(geodist,phi_mat_lumD,[],colormat,'filled','b');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Phase (rad)','Interpreter','latex');
% title(['Lumen intensity sum phase progression, v = ',num2str(round(resultsmat(1).phasevec_I(7),2)),' mm/s, $R^2$ = ',num2str(round(resultsmat(1).phasevec_I(2)^2,2))],'Interpreter','latex');
xplot = 0:0.001:max(geodist); yplot = resultsmat(1).phasevec_I(1)*xplot + resultsmat(1).phasevec_I(8);
% hold on; plot(xplot,yplot,'k');
%Calcium progression
% subplot(5,2,[7 8]);
scatter(geodist,phi_mat_Ca,[],colormat,'filled','g');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Phase (rad)','Interpreter','latex');
% title(['Calcium phase progression, v = ',num2str(round(resultsmat(1).phasevec_Ca(7),2)),' mm/s, $R^2$ = ',num2str(round(resultsmat(1).phasevec_Ca(2)^2,2))],'Interpreter','latex');
xplot = 0:0.001:max(geodist); yplot = resultsmat(1).phasevec_Ca(1)*xplot + resultsmat(1).phasevec_Ca(8);
% hold on; plot(xplot,yplot,'k');
legend({'D FWHM','$\Sigma$Intensity','Calcium'},'Interpreter','latex')

%Plot log-log slopes at each location
subplot(5,2,[7 8])
scatter(geodist,I_DSlope,'filled','k');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Slope log($\Sigma$I) vs log(D)','Interpreter','latex');
yline(2,'Alpha',0.3); ylim([0 4]);

%Plot Ca-Intensity lags
subplot(5,2,9);
scatter(geodist,I_Ca_lagvec,[],colormat,'filled');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Lag (seconds)','Interpreter','latex');
title('Calcium to Integrated Intensity lag','Interpreter','latex');
ylim([1 4]);
%Plot Ca-Diameter lags
subplot(5,2,10);
scatter(geodist,D_Ca_lagvec,[],colormat,'filled');
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Lag (seconds)','Interpreter','latex');
title('Calcium to FWHM Diameter lag','Interpreter','latex');
ylim([1 4]);

sgtitle(['JD221223F2 ',roiname],'Interpreter', 'none','FontSize',12);
print(fig,[roiname,'_PhaseProgression'],'-dpdf','-r0','-painters')

%% Plot intensity vs diameter (logscale)

fig=figure('units','inches','outerposition',[0 0 8.5 11]); hold on;
s1 = subplot(10,2,[1 2]);


plotlumD = 1;
for i = 1:size(diam,1)/2
% for i = size(diam,1)
    xfit = log10(diam_lp(i,:));
    yfit = log10(lumwave_lp(i,:));
    meany = mean(yfit);
%     xfit = log10(diam(i,:));
%     yfit = log10(lumwave(i,:));
    fittmp = polyfit(xfit,yfit,1);
    I_DSlope(i) = fittmp(1);
    I_DConst(i) = fittmp(2);
        y_est = polyval(fittmp,xfit);
    SSR = sum((yfit - y_est).^2);
    SST = sum((yfit - meany).^2);
    I_DR2(i) = 1-SSR/SST;
    if plotlumD == 1
        subplot(10,2,i);
        scatter(xfit,yfit,'filled');
        title(['Segment ',num2str(i),', slope = ',num2str(round(I_DSlope(i),2))],'Interpreter','latex');
        hold on; plot(xfit,y_est,'r--');
        ylim([min(min(log10(lumwave_lp(:,:)))),max(max(log10(lumwave_lp(:,:))))])
        xlim([min(min(log10(diam_lp(:,:)))),max(max(log10(diam_lp(:,:))))])
        ylabel('log10($\Sigma$I)','Interpreter','latex');
%         if i == 15 || i == 14
%         xlabel('log10($D_{FWHM}$)','Interpreter','latex');
%         end
        %     ylim([min(min(log10(lumwave(:,:)))),max(max(log10(lumwave(:,:))))])
        %     xlim([min(min(log10(diam(:,:)))),max(max(log10(diam(:,:))))])
    end
end

sgtitle(['JD221223F2 ',roiname,' Intensity vs FWHM Diameter'],'Interpreter', 'none','FontSize',12);
print(fig,[roiname,'_IntensityScaling'],'-dpdf','-r0','-vector')


