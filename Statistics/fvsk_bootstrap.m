%fvsk_bootstrap.m




function [b_bootstrap] = fvsk_bootstrap(f,k,wts)

trials = 10000; %Number of bootstrap trials
b_bootstrap = zeros(trials,1);
samples = numel(f);
data = [f,k,wts];
parfor i = 1:trials
    data_b = datasample(data,samples,1);
%     [wtdmedf,wtdmedk] =
%     weightedMedFit(data_b(:,2),data_b(:,1),data_b(:,3),10); %For arbitrary (continuous, rest state) peak vasomotor frequencies.
    [wtdmedf,wtdmedk] = weightedMedFit_discretefreq(data_b(:,2),data_b(:,1),data_b(:,3)); %For discrete (stim) frequencies
    [b,~,~] = FitSlope_NoIntercept(wtdmedk,wtdmedf);
    b_bootstrap(i) = b;
end


end