%Analyze All Vessels and Neurons
%Combine vasomotor phase gradient and frequency results across animals,
%trials, and vessels. Filter based on vessel length and presence of
%traveling waves. Make summary plots. This script is used after widefield_ExtractVesNeuPhaseGrad.m 

clear
close all
animal = 'TB013120M2';
cd(['AllSegments\',animal]); %Change directory to animal folder with extracted ves/neu phase gradient results.

%% Make and save combined phasevecs for each animal, without filtering for t value (can do that later if needed)
filesv = dir('**\*results.mat')
filesn = dir('**\*neural_results.mat')
todel = zeros(length(filesv),1);
for i=1:length(filesv) 
    todel(i) = contains(filesv(i).name,'neural'); 
end
filesv(logical(todel)) = [];

clearvars todel todeln
same = zeros(length(filesn),1);
for i=1:length(filesn)
    if strcmp(extractBefore(filesn(i).name,'neural_results.mat'),extractBefore(filesv(i).name,'results.mat'))
        same(i) = 1;
    end
end
filecheck = sum(same);
if filecheck == length(filesn)
    disp('V N files match')
    pause()
end
%%
for i=1:length(filesv)
    %Load vessel and neuron phase progression results
    data_folderv = filesv(i).folder;
    tmpv = load([data_folderv,'\',filesv(i).name]);
    link_results_struct = tmpv.link_results_struct;
    pix_mm = link_results_struct(1).phasevec(1,5)/link_results_struct(1).link_lengths_mm(1);
    data_foldern = filesn(i).folder;
    tmpn = load([data_foldern,'\',filesn(i).name]);
    link_results_struct_neu = tmpn.link_results_struct_neu;
    
    if i==1
        phasevec_combv = link_results_struct(1).phasevec;
        phasevec_combv = [phasevec_combv,link_results_struct(1).t_test_mat]; %Add t_test_mat to the end of phasevec
        phasevec_combn = link_results_struct_neu(1).n_phasevec;
        phasevec_combn = [phasevec_combn,link_results_struct_neu(1).neu_t_test_mat]; %Add t_test_mat to the end of phasevec
        
        [segs_link] = fun_calcsegs_link(link_results_struct);
        phasevec_combv(:,end+1) = segs_link; %Segs per link
        size_trial = size(phasevec_combv,1);
        phasevec_combv(1:size_trial,end+1) = pix_mm;

        [segs_linkn] = fun_calcsegs_linkn(link_results_struct_neu);
        phasevec_combn(:,end+1) = segs_linkn;
        size_trialn = size(phasevec_combn,1);
        phasevec_combn(1:size_trialn,end+1) = pix_mm;
        
        nantodel = isnan(phasevec_combv(:,1));
        t_todel = zeros(size(phasevec_combv(:,1),1),1);
        for j=1:length(t_todel)
            if abs(phasevec_combv(j,9))<phasevec_combv(j,10) %If |t|<tcrit can't reject H0
                t_todel(j) = 1;
            end
        end
        length_vec = phasevec_combv(:,5);
        length_vec = length_vec/pix_mm;
        mlength_mm = 0.75; %Pixels/mm
        todelv = nantodel; %length_todel %| t_todel; %Not deleting entries for length or t-value when making combined results matrix
        
        nantodeln = isnan(phasevec_combn(:,1));
        t_todeln = zeros(size(phasevec_combn(:,1),1),1);
        for j=1:length(t_todeln)
            if abs(phasevec_combn(j,9))<phasevec_combn(j,10) %If |t|<tcrit can't reject H0
                t_todeln(j) = 1;
            end
        end
        length_vecn = phasevec_combn(:,5);
        length_vecn = length_vecn/pix_mm;
        mlength_mm = 0.75; %Pixels/mm
%         length_todeln = length_vecn<mlength_mm;
        todeln = nantodeln;
        todel = or(todelv, todeln);
        
        phasevec_combv(todel,:) = []; %Delete NaN entries in phasevecs (links missing phase calculations)
        phasevec_combn(todel,:) = []; %Delete NaN entries in phasevecs (links missing phase calculations)
        phasevec_combv2 = phasevec_combv;
        phasevec_combn2 = phasevec_combn;
        
        i
    else
        phasevec_combv_tmp = link_results_struct(1).phasevec;
        size_trial = size(phasevec_combv_tmp,1);
        phasevec_combv_tmp = [phasevec_combv_tmp,link_results_struct(1).t_test_mat];
        phasevec_combn_tmp = link_results_struct_neu(1).n_phasevec;
        size_trialn = size(phasevec_combn_tmp,1);
        phasevec_combn_tmp = [phasevec_combn_tmp,link_results_struct_neu(1).neu_t_test_mat];
        
        [segs_link] = fun_calcsegs_link(link_results_struct);
        phasevec_combv_tmp = [phasevec_combv_tmp,segs_link];
        phasevec_combv_tmp(:,end+1) = pix_mm;
        phasevec_combv = [phasevec_combv;phasevec_combv_tmp]; %Concatenate with previous trials phasevecs.
        
        [segs_linkn] = fun_calcsegs_linkn(link_results_struct_neu);
        phasevec_combn_tmp(:,end+1) = segs_linkn;
        size_trialn = size(phasevec_combn,1);
        phasevec_combn_tmp(:,end+1) = pix_mm;
        phasevec_combn = [phasevec_combn;phasevec_combn_tmp]; %Concatenate with previous trials phasevecs.
        
        
        nantodel = isnan(phasevec_combv(:,1));
        t_todel = zeros(size(phasevec_combv(:,1),1),1);
        for j=1:length(t_todel)
            if abs(phasevec_combv(j,9))<phasevec_combv(j,10) %If t<tcrit can't reject H0
                t_todel(j) = 1;
            end
        end
        length_todel = zeros(size(phasevec_combv(:,1),1),1);
        length_vec = phasevec_combv(:,5);
        length_vec = length_vec/pix_mm;
        mlength_mm = 0.75; %Pixels/mm
        length_todel = length_vec<mlength_mm;
        todelv = nantodel; %| length_todel %| t_todel;
        
        nantodeln = isnan(phasevec_combn(:,1));
        t_todeln = zeros(size(phasevec_combn(:,1),1),1);
        for j=1:length(t_todeln)
            if abs(phasevec_combn(j,9))<phasevec_combn(j,10) %If t<tcrit can't reject H0
                t_todeln(j) = 1;
            end
        end
        length_todeln = zeros(size(phasevec_combn(:,1),1),1);
        length_vecn = phasevec_combn(:,5);
        length_vecn = length_vecn/pix_mm;
        mlength_mm = 0.75; %Pixels/mm
%         length_todeln = length_vecn < mlength_mm;
        todeln = nantodeln; % | length_todeln %| t_todeln;
        
        todel = or(todelv,todeln);
        
        phasevec_combv(todel,:) = []; %Delete NaN entries in phasevecs (links missing phase calculations)
        phasevec_combn(todel,:) = []; %Delete NaN entries in phasevecs (links missing phase calculations)
        phasevec_combv2 = phasevec_combv;
        phasevec_combn2 = phasevec_combn;

        i
    end
end
%% Save each animal's combined resutls
save(['BeforeStim',animal,'vesselfv_VesselcombPhaseVec.mat'],'phasevec_combv');
save(['BeforeStim',animal,'vesselfv_NeuroncombPhaseVec.mat'],'phasevec_combn');


%% Combine data from all animals
clear; clc; close all;

%Loop to load combined data
filesv = dir('**\*VesselcombPhaseVec_9_4.mat');
filesn = dir('**\*NeuroncombPhaseVec_9_4.mat');
todelv = zeros(length(filesv),1);
todeln = zeros(length(filesn),1);
for i=1:length(filesv) %This animal is open window, not thin-skull
    if contains(filesv(i).name,'JD221024F2')
        todelv(i) = 1;
    end
    if contains(filesn(i).name,'JD221024F2')
        todeln(i) = 1;
    end
end
filesv(logical(todelv)) = [];
filesn(logical(todeln)) = [];

for i=1:length(filesv)
    if i == 1
        pvtmp = load([filesv(i).folder,'\',filesv(i).name]);
        pvcomb = pvtmp.phasevec_combv;
        pntmp = load([filesn(i).folder,'\',filesn(i).name]);
        pvncomb = pntmp.phasevec_combn;
    else
        pvtmp = load([filesv(i).folder,'\',filesv(i).name]);
        pvcomb = [pvcomb;pvtmp.phasevec_combv];
        pntmp = load([filesn(i).folder,'\',filesn(i).name]);        
        pvncomb = [pvncomb;pntmp.phasevec_combn];
    end
    i
end

%% Filter the total dataset (all animals, all vessels) for length and phase vs. distance correlation.
%Filter for vessel length
minlength = 0.75;
lengthvec = pvcomb(:,5)./pvcomb(:,end);
length_todel = lengthvec < minlength;

%Re-calculate t mat for different alpha values
[t_test_mat] = fun_calc_ttestmat(pvcomb,0.01); %Calculate correlation coefficient t-value and tcrit for different significance levels
pvcomb(:,9:10) = t_test_mat;

t_todel = zeros(size(pvcomb,1),1);
for j=1:length(t_todel) %default t values calculated for alpha = 0.05
    if abs(pvcomb(j,9))<pvcomb(j,10)% || abs(pvncomb(j,9))<pvncomb(j,10) %If t<tcrit can't reject H0
        t_todel(j) = 1;
    end
end
t_todel = logical(t_todel);
todel = or(length_todel,t_todel); %Filter by vessel length and correlation between phase and distnace

pvcomb(todel,:) = [];
pvncomb(todel,:) = [];
% These are now saved for other analyses %

%% Plot
%F vs K (change alpha appropriately)
figure
alpha = 0.01;
scatter(abs(pvcomb(:,1)),pvcomb(:,3),'filled','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
xlim([0 2]); ylim([0 0.18]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex');
str = sprintf('Pial arteriole f vs k, %.0f vessels',length(pvcomb));
str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
str2 = sprintf('%.2f',alpha);
title({str,[str1,' $\alpha$ = ',str2]},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
print(gcf, '-depsc2', 'fvsk_75length_01Tfilt');
savefig('fvsk_75length_01Tfilt.fig');

%R2 vs k
figure
alpha = 0.01;
scatter(abs(pvcomb(:,1)),pvcomb(:,2).^2,'filled','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
xlim([0 1.5]); ylim([0 1]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Variance explained, R2','Interpreter','latex');
str = sprintf('Pial arteriole f vs k, 24 animals, %.0f vessels',length(pvcomb));
str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
str2 = sprintf('%.2f',alpha);
title({str,[str1,' $\alpha$ = ',str2]},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
print(gcf, '-depsc2', 'R2vsk_75length_01Tfilt');

%Fvaso histogram
figure
h1 = histogram(pvcomb(:,3),'BinWidth',0.01);
h1.Orientation = 'horizontal';
h1.FaceAlpha = 1;
ylim([0 0.18])
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex','FontSize',14);
xlabel('Count','Interpreter','latex','FontSize',14);
str = sprintf('%.0f vessels',length(pvcomb));
str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
str2 = sprintf('%.2f',alpha);
title({['Vasomotion Frequency Distribution ',str],[str1,' $\alpha$ = ',str2]},'Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11
print(gcf, '-depsc2', 'VasoFreq_MarginalHist_75length_01Tfilt');
savefig('VasoFreq_MarginalHist_75length_01Tfilt.fig')

[f,x] = ecdf(pvcomb(:,3));
figure
plot(f,x)
ylim([0 0.18])
ylabel('Frequency (Hz)', 'Interpreter','latex')
xlabel('Probability','Interpreter','latex')
title('Vasomotion Frequency Cumulative Distribution, Alpha = 0.01','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
print(gcf, '-depsc2', 'VasoFreq_MarginalCumulativeDist_Alpha0_01');
savefig('VasoFreq_MarginalCumulativeDist_Alpha0_01.fig')

%phase gradient histogram
figure
h1 = histogram(abs(pvcomb(:,1)),'BinWidth',0.05);
h1.FaceAlpha = 1;
xlim([0 2])
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex','FontSize',14);
ylabel('Count','Interpreter','latex','FontSize',14);
str = sprintf('%.0f vessels',length(pvcomb));
str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
str2 = sprintf('%.2f',alpha);
title({['Phase gradient Distribution ',str],[str1,' $\alpha$ = ',str2]},'Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11
print(gcf, '-depsc2', 'PhaseGrad_MarginalHist_75length_01Tfilt');
savefig('PhaseGrad_MarginalHist_75length_01Tfilt.fig')

figure
ecdf(abs(pvcomb(:,1)))
xlim([0 2])
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Probability','Interpreter','latex')
str = sprintf('%.0f vessels',length(pvcomb));
str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
str2 = sprintf('%.2f',alpha);
title({['Phase gradient Cumulative Distribution ',str],[str1,' $\alpha$ = ',str2]},'Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11
print(gcf, '-depsc2', 'PhaseGrad_MarginalCumulativeDist_Alpha0_01');
savefig('PhaseGrad_MarginalCumulativeDist_Alpha0_01.fig')

%Plot kv vs kn
%These should be filtered for both vessel and neuronal t-values
figure
scatter(pvncomb(:,1),pvcomb(:,1),18,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
x = linspace(-1.5,1.5,100);
y=x; hold on; plot(x,y,'Color',[1,0,0,0.1]);
xline(0); yline(0);
daspect([1,1,1]); grid on
xlabel('k neurons (rad/mm)','Interpreter','latex')
ylabel('k vessels (rad/mm)','Interpreter','latex')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
title('kv vs kn','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
print(gcf, '-depsc2', 'kvsk_75length');
savefig('kvsk_75length.fig')
