% preprocessdiam_VisCrosslines.m
% Detect outliers by manually inspecting each line's signal.

function [todel,diam_mf,diam_lp,diam_mean,diam_mean_mf,meand,filtparams] = preprocessdiam_VisCrosslines(diam,rate,inds,mask)
    %First identify cross sections that should be excluded
    dt = 1/rate;
    time = dt:dt:(size(diam,1)*dt);

    meand = mean(diam,1);
    meanreps = repmat(meand,[size(diam,1),1]);
    diam_mean = diam - meanreps;

    params.Fs = rate; %Interpolated rate is twice actual single depth rate
    params.pad = 2;
    params.fpass = [0 params.Fs/2]; %Hz, default is [0 Fs/2]
    params.err   = [2 .05];
    params.trialave = 0;
    T = time(end);
    BW = 0.02; %600s trial -> TBW =~ 7
    if T < 300
        BW = 0.04;
    end
    params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers
    addpath(genpath('chronux_2_12'))
    [S,f] = mtspectrumc(diam_mean,params);

    numpts = size(diam,2);
    colormat = zeros(numpts,3); colormat(:,1) = 0; colormat(:,3) = 1; colormat(:,2) = flip((1/numpts):(1/numpts):1);
    coloriter = 1;


    isusable = zeros(numpts,1); %Keep lines based on data quality (usually poor quality is due to underlying vessel)

    h1 = openfig('Ca_Lumen_imfust.fig','reuse');
    ax1 = gca;
    fig1 = get(ax1,'children');

    fig=figure('units','inches','outerposition',[0 0 8.5 11]); hold on;
    s1 = subplot(3,2,[1 2]);
    copyobj(fig1,s1);
    ylim([0 size(mask,1)]); xlim([0 size(mask,2)/4]);
    daspect([1,4,1]);
    axis off

    for i = 1:numpts
        axes(s1); hold on;

        linenum = i*10-5;
        ind1 = inds{1,linenum}(1);
        ind2 = inds{1,linenum}(end);
        [row1,col1] = ind2sub(size(mask),ind1);
        [row2,col2] = ind2sub(size(mask),ind2);
        line([col1,col2]/4,[row1,row2],'Color',colormat(coloriter,:)); 
        coloriter = coloriter + 1;


        subplot(3,2,[3 4]);
        plot(time,diam(:,i)); title(num2str(i)); ylim([min(diam(:)),max(diam(:))]); 
        subplot(3,2,[5 6]);
        plot(f,log10(S(:,i)/sum(S(:,i))),'b'); xlim([0 1]); ylim([min(log10(S(:,i)/sum(S(:,i)))) max(log10(mean(S,2)/sum(mean(S,2))))]);
        hold on; plot(f,log10(mean(S,2)/sum(mean(S,2))),'k');
        pause();
        isusable(i) = input('Is Point Usable? 1 Yes ');
        hold off
    end
    todel = ~isusable;
    
    diam(:,todel) = []; %Keep only points with no errors in detection
   
    %Median filter and output diameter data to use for time analysis 
    mfrange = 5;
    diam_mf = medfilt1(diam,5);
    diam_mean_mf = medfilt1(diam_mean,mfrange);

    FiltOrder = 300; SBFreq = 0.3; PBFreq = 0.25;
    lpFilt = designfilt('lowpassfir','FilterOrder',FiltOrder,'StopbandFrequency',SBFreq,'PassbandFrequency',PBFreq,'SampleRate',rate); %Bandpass Filter, requires large filter order, order must be 3x #points
    diam_lp = filtfilt(lpFilt,diam_mean);

    filtparams = struct();
    filtparams.mfrange = mfrange;
    filtparams.FiltOrder = FiltOrder;
    filtparams.SBFreq = SBFreq;
    filtparams.PBFreq = PBFreq;

end