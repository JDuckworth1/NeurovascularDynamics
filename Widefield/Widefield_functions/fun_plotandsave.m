% fun_saveandplot.m

function [] = fun_plotandsave(plotQ,link_results_struct,link_results_struct_neu,t_test_mat,neu_t_test_mat,vsl_graph,max_vec,n_max_vec,pix_mm,phasevec,n_phasevec,f,animal,wind)

if plotQ == 1
    figure
    h2 = histogram(max_vec/pix_mm,20);
    h2.BinWidth = 0.1;
    title('Vessel Length distribution after modifying graph')
    xlabel('Length (mm)')
    ylabel('Count')
    xlim([0 5])
    savefig('Vessel_LengthDist_ModifiedGraph.fig')

    figure
    h2 = histogram(n_max_vec/pix_mm,20);
    h2.BinWidth = 0.1;
    title('Neural Length distribution after modifying graph')
    xlabel('Length (mm)')
    ylabel('Count')
    xlim([0 5])
    savefig('Neural_LengthDist_ModifiedGraph.fig')

    figure
    scatter(max_vec/pix_mm,abs(t_test_mat(:,1))-t_test_mat(:,2),5,'filled')
    yline(0);xline(0.75);
    ylabel('Abs(t)-t_{crit}')
    xlabel('Vessel Segment Length (mm)')
    ylim([-10 10])
    savefig('Vessel_Ttest_vsLength.fig')

    figure
    scatter(max_vec/pix_mm,abs(neu_t_test_mat(:,1))-neu_t_test_mat(:,2),5,'filled')
    yline(0);xline(0.75);
    ylabel('Abs(t)-t_{crit}')
    xlabel('Vessel Segment Length (mm)')
    ylim([-10 10])
    savefig('Neural_Ttest_vsVesselLength.fig')
end

%Filter results to get rid of NaN, t-val fails, length < 0.75 mm
%Filter by vessel length, vessel t-val
R_vec = link_results_struct(1).phasevec(:,2);
k_vec = abs(link_results_struct(1).phasevec(:,1));
k_vec_sgn = link_results_struct(1).phasevec(:,1);
k_vec_sgn_fv = link_results_struct(1).phasevec_fv(:,1);
k_vec_sgn_fs = link_results_struct(1).phasevec_fs(:,1);
k_vec_sgn_fs2 = link_results_struct(1).phasevec_fs2(:,1);
length_vec = link_results_struct(1).link_lengths_mm';
length_vec = length_vec';
f_vec = link_results_struct(1).f_peak';
spec_vec = link_results_struct(1).AvgSpec;
phasevec2 = phasevec;
t_vec = t_test_mat;
cohvec = link_results_struct(1).C12_f;

n_R_vec = link_results_struct_neu(1).n_phasevec(:,2);
n_k_vec = abs(link_results_struct_neu(1).n_phasevec(:,1));
n_k_vec_sgn = link_results_struct_neu(1).n_phasevec(:,1);
n_k_vec_sgn_fv = link_results_struct_neu(1).n_phasevec_fv(:,1);
n_k_vec_sgn_fs = link_results_struct_neu(1).n_phasevec_fs(:,1);
n_k_vec_sgn_fs2 = link_results_struct_neu(1).n_phasevec_fs2(:,1);
n_length_vec = link_results_struct_neu(1).n_link_lengths_mm';
n_length_vec = n_length_vec';
n_phasevec2 = n_phasevec;
n_t_vec = neu_t_test_mat;

%Get rid of NaN entries
nantodel = isnan(k_vec);
n_nantodel = isnan(n_k_vec);
tot_nantodel = or(nantodel,n_nantodel); %To make sure we're deleting the same vessels from both vessel and neurons

k_vec(tot_nantodel) = [];
k_vec_sgn(tot_nantodel) = [];
k_vec_sgn_fv(tot_nantodel) = [];
k_vec_sgn_fs(tot_nantodel) = [];
k_vec_sgn_fs2(tot_nantodel) = [];
R_vec(tot_nantodel) = [];
length_vec(tot_nantodel) = [];
f_vec(tot_nantodel) = [];
spec_vec(:,tot_nantodel) = [];
phasevec2(tot_nantodel,:) = [];
t_vec(tot_nantodel,:) = [];
cohvec(:,tot_nantodel) = [];
rawnumsegs = length(k_vec);
tot_remain = [length(vsl_graph.link.cc_ind),length(k_vec)]

n_k_vec(tot_nantodel) = [];
n_k_vec_sgn(tot_nantodel) = [];
n_k_vec_sgn_fv(tot_nantodel) = [];
n_k_vec_sgn_fs(tot_nantodel) = [];
n_k_vec_sgn_fs2(tot_nantodel) = [];
n_R_vec(tot_nantodel) = [];
n_length_vec(tot_nantodel) = [];
n_phasevec2(tot_nantodel,:) = [];
n_t_vec(tot_nantodel,:) = [];

%Get rid of entries with vessel abs(t)<tcrit
t_todel = zeros(size(t_vec,1),1);
for i=1:length(t_todel)
    if abs(t_vec(i,1))<t_vec(i,2) %If t<tcrit
        t_todel(i) = 1;
    end
end
t_todel = logical(t_todel);
k_vec(t_todel) = [];
k_vec_sgn(t_todel) = [];
k_vec_sgn_fv(t_todel) = [];
k_vec_sgn_fs(t_todel) = [];
k_vec_sgn_fs2(t_todel) = [];
R_vec(t_todel) = [];
length_vec(t_todel) = [];
f_vec(t_todel) = [];
spec_vec(:,t_todel) = [];
phasevec2(t_todel,:) = [];
t_vec(t_todel,:) = [];
n_t_vec(t_todel,:) = [];
cohvec(:,t_todel) = [];
tot_remain = [length(vsl_graph.link.cc_ind),length(k_vec)]

n_k_vec(t_todel) = [];
n_k_vec_sgn(t_todel) = [];
n_k_vec_sgn_fv(t_todel) = [];
n_k_vec_sgn_fs(t_todel) = [];
n_k_vec_sgn_fs2(t_todel) = [];
n_R_vec(t_todel) = [];
n_length_vec(t_todel) = [];
n_phasevec2(t_todel,:) = [];

%Filter based on segment length
mlength_mm = 0.75; %mm
length_todel = length_vec<mlength_mm;
excluded_links = sum(length_todel) %Output how many were filtered out
k_vec(length_todel) = [];
k_vec_sgn(length_todel) = [];
k_vec_sgn_fv(length_todel) = [];
k_vec_sgn_fs(length_todel) = [];
k_vec_sgn_fs2(length_todel) = [];
R_vec(length_todel) = [];
length_vec(length_todel) = [];
f_vec(length_todel) = [];
spec_vec(:,length_todel) = [];
phasevec2(length_todel,:) = [];
n_t_vec(length_todel,:) = [];
cohvec(:,length_todel) = [];
numsegments = length(k_vec);
tot_remain = [length(vsl_graph.link.cc_ind),length(k_vec)]

n_k_vec(length_todel) = [];
n_k_vec_sgn(length_todel) = [];
n_k_vec_sgn_fv(length_todel) = [];
n_k_vec_sgn_fs(length_todel) = [];
n_k_vec_sgn_fs2(length_todel) = [];
n_R_vec(length_todel) = [];
n_length_vec(length_todel) = [];
n_phasevec2(length_todel,:) = [];
n_t_todel = abs(n_t_vec(:,1))<n_t_vec(:,2);

sig = phasevec2(:,4)*pix_mm;
numerator = sum((1./(sig.^2)).*abs(phasevec2(:,1)));
denom = sum(1./(sig.^2));
phasegrad_wtd = numerator/denom
% link_results_struct(1).phasegrad_wtd = phasegrad_wtd;
sig = n_phasevec2(:,4)*pix_mm; %Only ~15 entries after filtering above
numerator = sum((1./(sig.^2)).*abs(n_phasevec2(:,1)));
denom = sum(1./(sig.^2));
phasegrad_wtd = numerator/denom
% link_results_struct_neu(1).phasegrad_wtd = phasegrad_wtd;

if plotQ == 1
    %Plot results (total and filtered)
    figure %Plot phase gradient histogram
    h3 = histogram(k_vec);
    h3.BinWidth = 0.05;
    xlim([0 3])
    xlabel('Phase Gradient (rad/mm)')
    ylabel('Count')
    str = sprintf('Magnitude Vessel Phase Gradient Distribution %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title(str)
    savefig('VesselPhaseGradHist.fig')
    figure %Neurons
    h3 = histogram(n_k_vec);
    h3.BinWidth = 0.05;
    xlim([0 3])
    xlabel('Neuron Phase Gradient (rad/mm)')
    ylabel('Count')
    str = sprintf('Magnitude Neuron Phase Gradient Distribution %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title(str)
    savefig('NeuralPhaseGradHist.fig')

    figure
    histogram(f_vec,20)
    xlabel('Frequency (Hz)')
    ylabel('Count')
    xlim([0 0.2])
    savefig('VasofreqHist.fig')

    figure
    subplot(2,1,1)
    for i=1:size(spec_vec,2)
        plot(f,log10(spec_vec(:,i)))
        hold on
    end
    xlabel('Frequency (Hz)')
    ylabel('log10 Power')
    subplot(2,1,2)
    plot(f,log10(mean(spec_vec,2))); grid on;
    xlabel('Frequency (Hz)')
    ylabel('log10 Power')
    savefig('AvgSegSpectra.fig')

    figure %R^2 vs Phase Gradient, Vessels
    subplot(1,2,2)
    scatter(k_vec,R_vec.^2,5,'filled')
    ylabel('R^2')
    xlabel('Vessel Phase Gradient k (rad/mm)')
    str = sprintf('%d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title([animal,' abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 2])
    ylim([0 1])
    subplot(1,2,1)
    scatter(abs(link_results_struct(1).phasevec(:,1)),link_results_struct(1).phasevec(:,2).^2,5,'filled')
    ylabel('R^2')
    xlabel('Vessel Phase Gradient k (rad/mm)')
    str = sprintf('%d Vessel Segs, No Minimum Length',rawnumsegs);
    title([animal,' abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 2])
    ylim([0 1])
    savefig('VesselR2vsPhaseGrad.fig')
    figure %R^2 vs Phase Gradient, Neurons
    subplot(1,2,2)
    scatter(n_k_vec,n_R_vec.^2,5,'filled')
    ylabel('R^2')
    xlabel('Neuron Phase Gradient k (rad/mm)')
    str = sprintf('%d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title([animal,' abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 2])
    ylim([0 1])
    subplot(1,2,1)
    scatter(abs(link_results_struct_neu(1).n_phasevec(:,1)),link_results_struct_neu(1).n_phasevec(:,2).^2,5,'filled')
    ylabel('R^2')
    xlabel('Neuron Phase Gradient k (rad/mm)')
    str = sprintf('Neuons, %d Vessel Segs, No Minimum Length',rawnumsegs);
    title([animal,' abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 2])
    ylim([0 1])
    savefig('NeuralR2vsPhaseGrad.fig')

    figure %Phase Grad vs VasoFreq
    subplot(1,2,2)
    scatter(f_vec,k_vec,5,'filled')
    ylabel('Phase Gradient k (rad/mm)')
    xlabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('%d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title([animal,' abs($\phi$/$\l$) vs fv ',str],'Interpreter','latex')
    ylim([0 4])
    xlim(wind)
    subplot(1,2,1)
    scatter(link_results_struct(1).f_peak,abs(link_results_struct(1).phasevec(:,1)),5,'filled')
    ylabel('Phase Gradient k (rad/mm)')
    xlabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('%d Vessel Segs, No Minimum Length',rawnumsegs);
    title([animal,' abs($\phi$/$\l$) vs fv ',str],'Interpreter','latex')
    ylim([0 4])
    xlim(wind)
    savefig('VesselKvsVasoFreq.fig')
    figure %Phase Grad vs VasoFreq, Neurons
    subplot(1,2,2)
    scatter(f_vec,n_k_vec,5,'filled')
    ylabel('Neuron Phase Gradient k (rad/mm)')
    xlabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('Neurons, %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title([animal,' abs($\phi$/$\l$) vs fv ',str],'Interpreter','latex')
    ylim([0 4])
    xlim(wind)
    subplot(1,2,1)
    scatter(link_results_struct_neu(1).n_phasevec(:,3),abs(link_results_struct_neu(1).n_phasevec(:,1)),5,'filled')
    ylabel('Neuron Phase Gradient k (rad/mm)')
    xlabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('%d Vessel Segs, No Minimum Length',rawnumsegs);
    title([animal,' abs($\phi$/$\l$) vs fv ',str],'Interpreter','latex')
    ylim([0 4])
    xlim(wind)
    savefig('NeuralKvsVasoFreq.fig')

    figure %f vs k
    subplot(1,2,2)
    scatter(k_vec,f_vec,5,'filled')
    xlabel('Phase Gradient k (rad/mm)')
    ylabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('%d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title([animal,' fv vs abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 4])
    ylim(wind)
    subplot(1,2,1)
    scatter(abs(link_results_struct(1).phasevec(:,1)),link_results_struct(1).f_peak,5,'filled')
    xlabel('Phase Gradient k (rad/mm)')
    ylabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('%d Vessel Segs, No Minimum Length',rawnumsegs);
    title([animal,' fv vs abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 4])
    ylim(wind)
    savefig('VesselVasoFreqvsK.fig')
    figure %f vs k
    subplot(1,2,2)
    scatter(n_k_vec,f_vec,5,'filled')
    xlabel('Neural Phase Gradient k (rad/mm)')
    ylabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('Neurons, %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    title([animal,' fv vs abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 4])
    ylim(wind)
    subplot(1,2,1)
    scatter(abs(link_results_struct_neu(1).n_phasevec(:,1)),link_results_struct_neu(1).n_phasevec(:,3),5,'filled')
    xlabel('Neural Phase Gradient k (rad/mm)')
    ylabel('Segment Vasomotion Frequency (1/s)')
    str = sprintf('%d Vessel Segs, No Minimum Length',rawnumsegs);
    title([animal,' fv vs abs($\phi$/$\l$) ',str],'Interpreter','latex')
    xlim([0 4])
    ylim(wind)
    savefig('NeuralVasoFreqvsK.fig')

    %Vessel k vs Neural k
    % figure
    % scatter(n_k_vec_sgn,k_vec_sgn,'filled')
    % xlabel('k - neurons (rad/mm)','Interpreter','latex')
    % ylabel('k - vessels (rad/mm)','Interpreter','latex')
    % str = sprintf('vessel vs neuron phase gradient %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    % title(str,'Interpreter','latex')
    % grid on
    % xlim([-1 1])
    % ylim([-1 1])
    % xline(0);
    % yline(0);
    % hold on
    % n_k_vec(n_t_todel) = [];
    % n_k_vec_sgn(n_t_todel) = [];
    % k_vec_sgn(n_t_todel) = [];
    % scatter(n_k_vec_sgn,k_vec_sgn,'filled','r')
    % savefig('kvskScatter.fig')

    figure
    subplot(2,2,1)
    scatter(n_k_vec_sgn,k_vec_sgn,'filled')
    xlabel('k - neurons (rad/mm)','Interpreter','latex')
    ylabel('k - vessels (rad/mm)','Interpreter','latex')
    % str = sprintf('vessel vs neuron phase gradient %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    % title(str,'Interpreter','latex')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    xline(0);
    yline(0);
    hold on
    n_k_vec(n_t_todel) = [];
    n_k_vec_sgn(n_t_todel) = [];
    k_vec_sgn(n_t_todel) = [];
    % scatter(n_k_vec,k_vec_sgn,'filled','r') %Wasn't plotting signed version
    scatter(n_k_vec_sgn,k_vec_sgn,'filled','r')
    title('Phase at vessel Fv')
    subplot(2,2,2)
    scatter(n_k_vec_sgn_fv,k_vec_sgn_fv,'filled')
    xlabel('k - neurons (rad/mm)','Interpreter','latex')
    ylabel('k - vessels (rad/mm)','Interpreter','latex')
    % str = sprintf('vessel vs neuron phase gradient %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    % title(str,'Interpreter','latex')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    xline(0);
    yline(0);
    % hold on
    % n_k_vec(n_t_todel) = [];
    % n_k_vec_sgn_fv(n_t_todel) = [];
    % k_vec_sgn_fv(n_t_todel) = [];
    % scatter(n_k_vec,k_vec_sgn,'filled','r') %Wasn't plotting signed version
    % scatter(n_k_vec_sgn_fv,k_vec_sgn_fv,'filled','r')
    title('Phase at global Fv')
    subplot(2,2,3)
    scatter(n_k_vec_sgn_fs,k_vec_sgn_fs,'filled')
    xlabel('k - neurons (rad/mm)','Interpreter','latex')
    ylabel('k - vessels (rad/mm)','Interpreter','latex')
    % str = sprintf('vessel vs neuron phase gradient %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    % title(str,'Interpreter','latex')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    xline(0);
    yline(0);
    % hold on
    % n_k_vec(n_t_todel) = [];
    % n_k_vec_sgn_fs(n_t_todel) = [];
    % k_vec_sgn_fs(n_t_todel) = [];
    % scatter(n_k_vec,k_vec_sgn,'filled','r') %Wasn't plotting signed version
    % scatter(n_k_vec_sgn_fs,k_vec_sgn_fs,'filled','r')
    title('Phase at Fs')
    subplot(2,2,4)
    scatter(n_k_vec_sgn_fs2,k_vec_sgn_fs2,'filled')
    xlabel('k - neurons (rad/mm)','Interpreter','latex')
    ylabel('k - vessels (rad/mm)','Interpreter','latex')
    % str = sprintf('vessel vs neuron phase gradient %d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    % title(str,'Interpreter','latex')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    xline(0);
    yline(0);
    % hold on
    % n_k_vec(n_t_todel) = [];
    % n_k_vec_sgn_fs2(n_t_todel) = [];
    % k_vec_sgn_fs2(n_t_todel) = [];
    % scatter(n_k_vec,k_vec_sgn,'filled','r') %Wasn't plotting signed version
    % scatter(n_k_vec_sgn_fs2,k_vec_sgn_fs2,'filled','r')
    title('Phase at Fs 2nd Harmonic')
    savefig('kvskScatter.fig')

    % cohsdm = 1/sqrt(length(segs_link))*std(cohvec,0,2,'omitnan');
    cohsd = std(cohvec,0,2,'omitnan');
    figure
    plot(link_results_struct(1).freqforcoh,mean(abs(cohvec),2,'omitnan'));
    xlabel('Frequency (Hz)','Interpreter','latex')
    ylabel('Magnitude Coherence','Interpreter','latex')
    title('Average over vessels, magnitude coherence +/- SD','Interpreter','latex')
    hold on
    lh = plot(link_results_struct(1).freqforcoh,mean(abs(cohvec)+cohsd,2,'omitnan'),'r');
    lh2 = plot(link_results_struct(1).freqforcoh,mean(abs(cohvec)-cohsd,2,'omitnan'),'r');
    lh.Color = [0,0,0,0.3];
    lh2.Color = [0,0,0,0.3];
    savefig('AverageVNCoherenceMag.fig')
    pause(3);


    % figure
    % subplot(1,2,2)
    % scatter(f_vec,R_vec.^2,5,'filled')
    % ylabel('R^2')
    % xlabel('Segment Vasomotion Frequency (1/s)')
    % str = sprintf('%d Vessel Segs, %.2f mm Min Seg Length',numsegments,mlength_mm);
    % title([animal,' abs($\phi$/$\l$) vs fv ',str],'Interpreter','latex')
    % ylim([0 1])
    % xlim(wind)
    % subplot(1,2,1)
    % scatter(link_results_struct(1).f_peak,link_results_struct(1).phasevec(:,2).^2,5,'filled')
    % ylabel('R^2')
    % xlabel('Segment Vasomotion Frequency (1/s)')
    % str = sprintf('%d Vessel Segs, No Minimum Length',rawnumsegs);
    % title([animal,' abs($\phi$/$\l$) vs fv ',str],'Interpreter','latex')
    % ylim([0 1])
    % xlim(wind)
end









end