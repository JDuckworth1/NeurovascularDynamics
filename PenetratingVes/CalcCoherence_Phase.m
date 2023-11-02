% CalcCoherence_Phase.m
% Load extracted diameter at both depths. Perform outlier detection.
% Linearly interpolate each signal and calculate the coherence and phase
% between depths.

% Import data from each trial
clear; clc; close all;

trialnum = 1; %This is the experiment session/day
%Base folder varies, needs to be updated by user
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

% Set save location
save_folder = data_folder;

%%
stim = 'rest';
for file=1:length(files)

    [PA,depth1,depth2,~,pix_um] = fun_get_PAMag_depth(files(file).name);
    rate = 7.25; %Hz, constant for all trials.

    cd(save_folder)
    %Load struct with extracted diameters
    resultfiles = dir('*allstats.mat'); 
    tokeep = zeros(length(resultfiles),1);
    for i=1:length(resultfiles)
        if contains(resultfiles(i).name,['_PA',PA,'_'])
            tokeep(i) = 1;
        end
    end
    resultfiles(~logical(tokeep)) = [];

    isrest = zeros(length(resultfiles),1);
    isstim = zeros(length(resultfiles),1);
    for i=1:length(resultfiles)
        if contains(resultfiles(i).name,'stim') %Stim periods are differentiated in their struct's name.
            isstim(i) = 1;
        end
    end

    if strcmp(stim,'stim');
        resultfiles(~logical(isstim)) = [];
    end
    disp(resultfiles(1).name);
    disp(resultfiles(2).name);

    clearvars -except data_folder animal file files namestr depth_str depth1 depth2 PA stim pix_um rate resultfiles resultfilestmp analysisfolder phi_mat phi_struct oldfile oldpath oldphi_struct save_folder

    stats1 = load(resultfiles(1).name);
    stats2 = load(resultfiles(2).name);
    stats1 = stats1.allstats;
    stats2 = stats2.allstats;
    if str2double(stats1.depth) > str2double(stats2.depth) %Ensure stats1 is the shallow time series.
        stats1 = load(resultfiles(2).name);
        stats2 = load(resultfiles(1).name);
        stats1 = stats1.allstats;
        stats2 = stats2.allstats;
    end

   %% Perform outlier detection
   [stats1] = fun_outlier_rej(stats1,rate);
   [stats2] = fun_outlier_rej(stats2,rate);

    %% Calculate coherence and phase
    waves = stats1.RadondEq_Outl;
    waves = waves - mean(waves);
    waved = stats2.RadondEq_Outl;
    waved = waved - mean(waved);
    times = stats1.time;
    timed = stats2.time;
    if sum(times - timed) == 0 %Make sure that the time vectors reflect average frame times
        timed = times + (1/rate)/2;
        disp('Check times'); pause();
    end

    %Interpolate each signal to get equal time points
    tqs = times(1):(1/(2*rate)):times(end);
    tqd = timed(1):(1/(2*rate)):timed(end);
    waves_q_tmp = interp1(times,waves,tqs);
    waved_q_tmp = interp1(timed,waved,tqd);

    %Delete first wave1 value so that both start at t=0.75.
    waves_q = waves_q_tmp(2:end);
    tqs = tqs(2:end);
    %Delete last wave1 value to make both time series the same duration
    waved_q = waved_q_tmp;
    waved_q(end) = [];
    tqd(end) = [];

    params.Fs = rate*2; %Interpolated rate is twice actual single depth rate
    params.pad = 2;
    params.fpass = [0 params.Fs/4]; %Hz, default is [0 Fs/2]
    params.err   = [2 .05];
    params.trialave = 0;
    T = stats1.time(end);
    BW = 0.02; %600s trial -> TBW =~ 12
    params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers
    addpath(genpath('C:\chronux_2_12'))

    [S,f] = mtspectrumc(waves_q,params);
    if strcmp(stim,'rest')
        figure
        plot(f,log10(S)); xlim([0 1]);
        [wind,~] = ginput(2); %Manually define peak limits (exclude low freq shoulder)
        findf1 = round(wind(1),3);
        findf2 = round(wind(2),3);
        f1 = find(round(f,3)==round(findf1,3),1,'first'); %Hz lower bound
        f2 = find(round(f,3)==round(findf2,3),1,'first'); %Hz upper bound
        rmpath(genpath('C:\chronux_2_12'))
        [maxes,f_peaks] = findpeaks(S(f1:f2),f(f1:f2));
        maxloc = find(maxes==max(maxes));
        f_peak = f_peaks(maxloc)
        freq_findx = max(find(round(f,3)==round(f_peak,3)));
        addpath(genpath('C:\chronux_2_12'))
    else
        f_peak = 0.1; %Hz
        freq_findx = max(find(round(f,3)==round(f_peak,3)));
        addpath(genpath('C:\chronux_2_12'))
    end
    %Calcuate phase between deep and shallow, using interpolated data
    params_err = params;
    params_err.err = [2,0.05];
    [C,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(waved_q,waves_q,params_err);
    dof = 2*round(2*T*BW-1); %Single value dof (2*tapers)
    df = 1./((dof/2)-1);
    p = 0.05; %Confidence level
    confC = sqrt(1 - p.^df);

    %Testing different (constant) frequencies 
    findx_05 = max(find(round(f,4)==0.05));
    findx_1 =  max(find(round(f,4)==0.1));
    findx_15 = max(find(round(f,4)==0.15));

    figure
    plot(f,phi,'b'); hold on; plot(f,phi + 2*phistd','Color',[0,0,1,0.2]); plot(f,phi - 2*phistd','Color',[0,0,1,0.2]); xlim([0 0.5]); ylim([-pi pi]);
    yline(0,'Alpha',0.9); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('Relative Phase (rad)','Interpreter','latex');
    xline(f_peak,'-.','$f_v$','Interpreter','latex','Alpha',0.9,'Color','k');
    xline(0.05,'-.','$f_v$','Interpreter','latex','Alpha',0.5,'Color','r');
    xline(0.1,'-.','$f_v$','Interpreter','latex','Alpha',0.5,'Color','r');
    xline(0.15,'-.','$f_v$','Interpreter','latex','Alpha',0.5,'Color','r');

    %Calculate phase gradient
    dz = (str2num(depth2) - str2num(depth1))/1000; %mm depth difference (mm)
    dx = sqrt((stats1.vescentery - stats2.vescentery)^2 + (stats1.vescenterx - stats2.vescenterx)^2); %pixels, horizontal distance between vessel segments
    dx = dx/pix_um/1000; %mm
    ds = sqrt(dz^2 + dx^2); %mm, Lower estimate of vessel distance between two planes (not completely straight so still an underestimate)

    phi_manualf = phi(freq_findx);
    k_manualf = phi_manualf / ds; %rad/mm
    v_manualf = (2*pi*f_peak)/k_manualf; %velocity mm/s, assuming linear dispersion relation
    phi_05 = phi(findx_05);
    k_05 = phi_05 / ds; %rad/mm
    v_05 = (2*pi*0.05)/k_05; %velocity mm/s, assuming no dispersion
    phi_1 = phi(findx_1);
    k_1 = phi_1 / ds; %rad/mm
    v_1 = (2*pi*0.1)/k_1; %velocity mm/s, assuming no dispersion
    phi_15 = phi(findx_15);
    k_15 = phi_15 / ds; %rad/mm
    v_15 = (2*pi*0.15)/k_15; %velocity mm/s, assuming no dispersion

    str1 = sprintf('%.2f',phi_manualf);
    str2 = sprintf('%.2f',abs(k_manualf));
    str3 = sprintf('%.2f',v_manualf);
    title({'Relative Phase and 95 Percent CI',['$\phi_{fv}$ = ',str1,', $\mid k \mid$ = ',str2,' $v$ = ',str3]},'Interpreter','latex')

    if file == 1
        phi_struct = struct();
    end
    phi_struct(file).PA = PA;
    phi_vec_tmp = [phi_manualf;phi_05;phi_1;phi_15];
    k_vec_tmp = [k_manualf;k_05;k_1;k_15];
    phi_uncert_tmp = [phistd(freq_findx);phistd(findx_05);phistd(findx_1);phistd(findx_15)];
    fv_tmp = [f_peak;0.05;0.1;0.15];

    phi_mat = [phi_vec_tmp,k_vec_tmp,phi_uncert_tmp,fv_tmp];
    
    % For these trials, vessel location relative to the first PA were
    % recorded. For the others, locations were calculated using a
    % wide-field image and 2P frame scan for distance calibration. 
    if contains(animal,'PA8') || contains(animal,'PA9') || contains(animal,'PA10') || contains(animal,'PA11')
        disp(['PA ',PA])
        prompt = "Enter recorded x location ";
        phi_struct(file).xloc = input(prompt); %um relative to PA1
        prompt = "Enter recorded y location ";
        phi_struct(file).yloc = input(prompt); %um relative to PA1
    else %Distances have been calculated previously, load these results and assign them appropriately here.
        if file == 1
            disp('Choose old phi_struct')
            disp(animal)
            [oldfile,oldpath] = uigetfile
        end
       
        oldphi_mat = load([oldpath,'\',oldfile]) %Here we already found distnces using the widefield image
        oldphi_struct = oldphi_mat.phi_struct;
        todel = zeros(length(oldphi_struct),1);
        for j = 1:length(oldphi_struct)
            if ~strcmp(PA,oldphi_struct(j).PA)
                todel(j) = 1;
            end
        end
        oldphi_struct(logical(todel)) = [];
        phi_struct(file).xloc = oldphi_struct(1).xloc;
        phi_struct(file).yloc = oldphi_struct(1).yloc;
    end
    phi_struct(file).phi_mat = phi_mat; %phase, |k|, and uncertainties at peak vasomotor frequency. 
    phi_struct(file).f = f; %Peak vasomotor frequency
    phi_struct(file).C = C;
    phi_struct(file).Cerr = Cerr;
    phi_struct(file).confC = confC;

    phi_struct(file).Ss_Sd = S2(freq_findx)/S1(freq_findx);
    phi_struct(file).Ss = S2(freq_findx);
    phi_struct(file).Sd = S1(freq_findx);
    phi_struct(file).confC = confC;

end

cd(save_folder);
save(['phi_struct_stim_',animal,'.mat'],'phi_struct');


%% Make CorrMat for pair-wise correlation in vasomotor wave travel direction.
cd(save_folder);

numpairs = nchoosek(length(phi_struct),2)
CorrMat = zeros(numpairs,12);
counter = 1;
for i=1:length(phi_struct)-1
    v1 = 2*pi/(phi_struct(i).phi_mat(1,2))*phi_struct(i).phi_mat(1,4); %Using v = 2*pi*f/k, negative bottom to top, positive top to bottom
    Pnum1 = str2num(phi_struct(i).PA);
    Pos1x = phi_struct(i).xloc;
    Pos1y = phi_struct(i).yloc;
    phi1 = phi_struct(i).phi_mat(1,1); %phi
    dphi1 = phi_struct(i).phi_mat(1,3); %std (phi)

    freq_findx = max(find(round(phi_struct(i).f,3)==round(phi_struct(i).phi_mat(1,4),3)));
    Cerrfv1 = phi_struct(i).C(freq_findx); %Get coherence value at vasomotion frequency
    if isfield(phi_struct,'confC')
        confC1 = phi_struct(i).confC;
    else
        confC1 = 0.3568; %Same time length, tapers, confidence level used always 
    end

    for j=(i+1):length(phi_struct)
        v2 = 2*pi/(phi_struct(j).phi_mat(1,2))*phi_struct(j).phi_mat(1,4);
        Pnum2 = str2num(phi_struct(j).PA);
        Pos2x = phi_struct(j).xloc;
        Pos2y = phi_struct(j).yloc;

        CorrMat(counter,1) = Pnum1;
        CorrMat(counter,2) = Pnum2;
        CorrMat(counter,3) = v1;
        CorrMat(counter,4) = v2;

        if contains(animal,'PA8') || contains(animal,'PA9') || contains(animal,'PA10') || contains(animal,'PA11')
            CorrMat(counter,5) = sqrt((Pos2x-Pos1x)^2 + (Pos2y-Pos1y)^2) / 1000; %Pair Wise distance (mm) Use this for WT8-12 where locations are recorded
        else
            cd(analysisfolder)
            load('wdf_pix_mm_struct.mat');
            Wdf_pix_mm = wdf_pix_mm_struct.value;
            CorrMat(counter,5) = sqrt((Pos2x-Pos1x)^2 + (Pos2y-Pos1y)^2) / Wdf_pix_mm; %Pair Wise distance (mm)
        end

        CorrMat(counter,6) = phi1; %phase 1
        CorrMat(counter,7) = dphi1; %phase 1 standard dev
        CorrMat(counter,8) = phi_struct(j).phi_mat(1,1); %phase 2
        CorrMat(counter,9) = phi_struct(j).phi_mat(1,3); %phase 2 standard dev

        CorrMat(counter,10) = phi_struct(j).Ss; %Shallow power at vasomotor frequency
        CorrMat(counter,11) = phi_struct(j).Sd; %Deep power at vasomotor frequency

        freq_findx = max(find(round(phi_struct(j).f,3)==round(phi_struct(j).phi_mat(1,4),3)));
        Cerrfv2 = phi_struct(j).C(freq_findx);

        if isfield(phi_struct,'confC')
            confC2 = phi_struct(i).confC;
        else
            confC2 = 0.3568; %Same time length, tapers, confidence level used always
        end

        if Cerrfv1 > confC1 && Cerrfv2 > confC2
            CorrMat(counter,12) = 1; %We have significant coherence
        else
            CorrMat(counter,12) = 0; %We don't have significant coherence
        end

        counter = counter + 1;
    end

    KFmat(i,1) = phi_struct(i).phi_mat(1,4); %fv for this vessel (from top cross section)
    KFmat(i,2) = phi_struct(i).phi_mat(1,2); %k_fv for this vessel
    KFmat(i,3) = phi_struct(i).phi_mat(1,3)/abs(phi_struct(i).phi_mat(1,1))*abs(phi_struct(i).phi_mat(1,2)); %std k_fv (from std phi, assume no uncertainty in distance)
    KFmat(i,4) = phi_struct(i).Ss;
    KFmat(i,5) = phi_struct(i).Sd;

    %Get coherence at vasofreq (in phi_mat(1,4))
    freq_findx = max(find(round(phi_struct(i).f,3)==round(phi_struct(i).phi_mat(1,4),3)));
    Cerrfv = phi_struct(i).C(freq_findx);% Coherence value
    if isfield(phi_struct,'confC')
        confC = phi_struct(i).confC;
    else
        confC = 0.3568; %Same time length, tapers, confidence level used always 
    end
    if Cerrfv > confC
        KFmat(i,6) = 1; %We have significant coherence
    else
        KFmat(i,6) = 0; %We don't have significant coherence
    end

    if i==length(phi_struct)-1
        KFmat(i+1,1) = phi_struct(i+1).phi_mat(1,4); %fv for this vessel (from top cross section)
        KFmat(i+1,2) = phi_struct(i+1).phi_mat(1,2); %k_fv for this vessel
        KFmat(i+1,3) = phi_struct(i+1).phi_mat(1,3)/abs(phi_struct(i+1).phi_mat(1,1))*abs(phi_struct(i+1).phi_mat(1,2)); %std k_fv (from std phi, assume no uncertainty in distance)
        KFmat(i+1,4) = phi_struct(i+1).Ss;
        KFmat(i+1,5) = phi_struct(i+1).Sd;
        freq_findx = max(find(round(phi_struct(i+1).f,3)==round(phi_struct(i+1).phi_mat(1,4),3)));
        Cerrfv = phi_struct(i+1).C(freq_findx); %Test for mean coherence > Confidence level.

        if isfield(phi_struct,'confC')
            confC = phi_struct(i+1).confC;
        else
            confC = 0.3568; %Same time length, tapers, confidence level used always
        end
        if Cerrfv > confC
            KFmat(i+1,6) = 1; %We have significant coherence
        else
            KFmat(i+1,6) = 0; %We don't have significant coherence
        end
    end
end

cd(save_folder);
Corr_KF_struct(1).phi_struct = phi_struct;
Corr_KF_struct(1).CorrMat = CorrMat;
Corr_KF_struct(1).KFmat = KFmat;
save([animal,'_Corr_KF_struct.mat'],'Corr_KF_struct')



