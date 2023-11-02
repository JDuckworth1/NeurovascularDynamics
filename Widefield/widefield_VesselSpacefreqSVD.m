%Calculate the space-frequency svd for vessel data and save in toplot mat

clc;clear;close all;
%Change directories to the folder containing data
% files = dir('10.Oct.2023*\*meters.mat') %Get files by date (imaging session)
files = dir('*TB013120M4*\*meters.mat') %Get files by animal
str = 'v' %Compute for vessels or neurons? ('v'/'n')

%% Process Vascular Data
for file = 1:length(files)

%% Data loading
clearvars -except incH5 str nfiles file files batch upn unprocessed endfile startfile
cd(files(file).folder);
fname = files(file).name; 

filename = extractBefore(files(file).name,'_Parameters.mat')
data_folder = files(file).folder;
deter = 'p'
if deter == 'p'
    toplottmp = load(append(data_folder,'\',filename,'svd_puff.mat'));
elseif deter == 'b'
    toplottmp = load(append(data_folder,'\',filename,'svd_beforepuff.mat'));
end

toplot = toplottmp.toplot;
fpeak1 = toplot.f_peak(1);
if deter == 'p' || deter == 'a'
    btoplot = load((append(data_folder,'\',filename,'svd_beforepuff.mat')),'toplot');
    btoplot = btoplot.toplot;
    toplot.f_peak(1) = btoplot.f_peak(1);
    toplot.wind = btoplot.wind;
end

toplot.fname=files(file).name;
str1=["dualCH","PwrG80R","RedGreen","PowerG80R","dualLED"];
str2 =["G80_","singleCH","_Single_","5Power","SingleCH","G05"];
[NrOCh] = fun_get_NrOCh(str,fname,str1,str2); %Determine number of channels for loading

rmpath(genpath('C:\chronux_2_12'));
SensorData = toplottmp.SensorData;
Sensorprams = toplottmp.Sensorprams;
recprams = toplottmp.recprams;

[~,idx] = findpeaks(diff(SensorData.p95b),'MinPeakHeight',1.5);
% [~,pledidx] = findpeaks(diff(SensorData.stimLED),'MinPeakHeight',1);
[~,pledidx] = findpeaks(diff(SensorData.leftp),'MinPeakHeight',1);
pridx = pledidx;
% [~,stimidx] = findpeaks(diff(SensorData.stimLED),'MinPeakHeight',1);
[~,stimidx] = findpeaks(diff(SensorData.leftp),'MinPeakHeight',1);
aprams.rpuffdiff=pridx;
figure
hold on
if isfield(SensorData,'stimLED')==1;
    plot(SensorData.stimLED)
end

%------------------------------------------------------------
% % Load h5 file
%% Load wave files
if deter == 'p'
    wave=h5read([erase(files(file).name,'_Parameters.mat'),'svd_puff_wave.h5'],'/wave');
elseif deter == 'b'
    wave=h5read([erase(files(file).name,'_Parameters.mat'),'svd_beforepuff_wave.h5'],'/wave');
end

%% find vasomotion peak if havent already
if toplot.f_peak(1) == 0
        figure
        plot(log10(toplot.mpowr), 'k');
        prompt='Define Start and end of the peak search: ';
        [toplot.wind,~] = ginput(2);
        toplot.wind = fix(toplot.wind);
        [~, fpeak1] = max(toplot.mpowr(:, toplot.wind(1):toplot.wind(2)), [], 2);
        fpeak1 = toplot.f(fpeak1 + toplot.wind(1) - 1);
        toplot.f_peak(1)=fpeak1;
end
%% Parameters
%parameters
% Perform FFT
% Update half-bandwidth according to the rounded p value

[U,S,V]=svd(wave,0);
sig_modes=round(length(V)/2);
Un=single(U(:,1:sig_modes));
Sn=single(S(1:sig_modes,1:sig_modes));
Vn=single(V(:,1:sig_modes));
if str == 'n'
    clear U S V
end
clear wave

ext=extractBetween(files(file).name,'_0.','Hz');
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
padding_ratio = 2; %padded time base is 4 times higher then 2^nextpow2(size(wave,2));
toplot.Delta_f = 0.02; % Hz. denoted as half-bandwidth 
if toplot.f_peak(2)<0.05;
    toplot.Delta_f = 0.01
elseif size(toplot.f_peak,2)>3 && abs(toplot.f_peak(2)-toplot.f_peak(3))<2*toplot.Delta_f
    toplot.Delta_f=(abs(toplot.f_peak(2)-toplot.f_peak(3))-0.002)/2
end
num_pixel = size(Un,1);
num_frame= size(Vn,1);
padding_ratio = 2; %padded time base is 4 times higher then 2^nextpow2(size(wave,2));
num_frame_pad = (2 ^ ceil(log2(num_frame))) * padding_ratio;
toplot.pad = num_frame_pad;
p = round(num_frame * toplot.Delta_f / toplot.rate); % (num_frame / acq_frequency = acquisition time, i.e. T in the paper)
toplot.Delta_f = p * toplot.rate / num_frame;
disp(['Bandwidth = ', num2str(toplot.Delta_f), ' Hz'])
toplot.num_tapers = 2 * p - 1;
[slep,~] = dpss(num_frame, p, toplot.num_tapers);

%% Perform Coherence FFT
ext=extractBetween(files(file).name,'_0.','Hz');
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

% Update half-bandwidth according to the rounded p value
if str == 'v'
    if isunix
        addpath(genpath())
    else
        addpath(genpath('C:\chronux_2_12'))
    end
    ntapers = [(toplot.num_tapers+1)/2,toplot.num_tapers];
    Fs = toplot.rate;
    nfft = toplot.pad;

    tapers = dpsschk(ntapers,size(V,1),Fs);
    [f,findx] = getfgrid(Fs,nfft,[0,Fs/2]);

    %FFT
    clear taperedFFT S_tot wave z
    taperedFFT = complex(zeros(length(f),ntapers(2),size(V,2)));
    for i = 1:size(V,2)
        J = mtfftc(V(:,i),tapers,nfft,Fs);
        J = J(findx,:,:);
        taperedFFT(:,:,i) = J;
    end

    scoresC = U*S;
    clear U S i J tapers
    if isunix
        rmpath(genpath(''));
    else
        rmpath(genpath('C:\chronux_2_12'));
    end
end
%% Calculate Coherence versus frequency
plot_num_f_points = fix(4*toplot.rate/2 / toplot.Delta_f); % 4X oversampling
toplot.coherence = zeros(plot_num_f_points,1);
toplot.f_vector = linspace(0,toplot.rate/2,plot_num_f_points);
tic
interp_FFT = interp1(f,reshape(taperedFFT,size(taperedFFT,1),[]),toplot.f_vector);
interp_FFT = reshape(interp_FFT,[],size(taperedFFT,2),size(taperedFFT,3));
parfor i = 1:length(toplot.f_vector)
    if str == 'v'
        m = scoresC*squeeze(interp_FFT(i,:,:))';
    else
        m = scores*squeeze(interp_FFT(i,:,:))';
    end
    s = svd(m,0);
    coherence(i) = squeeze(s(1))^2./sum(s.^2); %Global coherence at each frequency in Observed brain dynamics p210
end
toplot.coherence = coherence;
toc
disp('Done')

clear i m s interp_FFT coherence
%% Make U and save Prameters
tic
close all
f_global = toplot.f_peak;
interp_FFT = interp1(f,reshape(taperedFFT,size(taperedFFT,1),[]),f_global);
interp_FFT = reshape(interp_FFT,[],size(taperedFFT,2),size(taperedFFT,3));
toplot.U = zeros(length(toplot.skel_label),toplot.num_tapers,length(toplot.f_peak));
for i = 1:length(toplot.f_peak)
    if str == 'v'
        m = scoresC*squeeze(interp_FFT(i,:,:))';
    else 
        m = scores*squeeze(interp_FFT(i,:,:))';
    end
    [u,~,~] = svd(m,0);
    toplot.U(:,:,i) = u(toplot.skel_label, :);
end
toplot.U=permute(toplot.U, [3,1,2]);
clear i m s interp_FFT f_global
toc

%% Save Data
if exist('SensorData','var')
    if str == 'v'
        if deter == 'p' %Make new .mat file
            save([erase(files(file).name,'_Parameters.mat'),'svd_puff.mat'], 'toplot' , 'recprams' ,'Sensorprams','SensorData', 'aprams', 'Vn','scores');
        elseif deter == 'b'
            save([erase(files(file).name,'_Parameters.mat'),'svd_beforepuff.mat'], 'toplot' , 'recprams' ,'Sensorprams','SensorData', 'aprams', 'Vn','scores');
        end
    elseif str =='n'
        save(files(file).name, 'toplot' , 'recprams' , 'Sensorprams','SensorData', 'aprams' ,'Vn','scores', 'rim'); %'-append');
    end
else
    save(files(file).name, 'toplot' , 'recprams' , 'Vn','scores');
end
disp('Done Saving'),
close all
%save([toplot.fname,'.mat'], 'toplot' , 'recprams' , 'Sensorprams','SensorData', 'aprams', 'rim');
end
disp('Vessels Done')
 