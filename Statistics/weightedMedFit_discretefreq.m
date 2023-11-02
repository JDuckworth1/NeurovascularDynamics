

function [fpts,wtdmedk] = weightedMedFit_discretefreq(kvec,fvec,wts)

% For stim data, we have a set of discrete frquencies. Binning them by
% cumulative percent causes problems because we'd have to split up points
% of the same frequency into different bins. So here we assign a median k
% for each frequency in fvec.

% fbins = binvec';
fpts = unique(fvec);
meank = zeros(size(fpts,1),1);
wtdmedk = zeros(size(fpts,1),1);
sek = zeros(size(fpts,1),1);
medk = zeros(size(fpts,1),1);
fplot = zeros(size(fpts,1),1);
num_bin = zeros(size(fpts,1),1);


for i=1:size(fpts,1)

    inbin = fvec == fpts(i);
    fvectmp = fvec(inbin);
    kvectmp = kvec(inbin);
    wtstmp = wts(inbin);
 
    meank(i) = sum(wtstmp.*kvectmp)/sum(wtstmp); %Weighted mean k within bin

    wtstot = sum(wtstmp);
    [kvectmpsort,sortinds] = sort(kvectmp); %Sort k in increasing order
    wtstmpsort1 = wtstmp(sortinds); %Sort weights by k order

    %To be the weighted-median point, sum of weights before and sum of
    %weights after must be <= half the sum(weights)
    cond1 = zeros(length(kvectmp),1);
    cond2 = zeros(length(kvectmp),1);
    for j = 1:length(kvectmp)
        if j ~= 1 && j ~= length(kvectmp)
            cond1(j) = sum(wtstmpsort1(1:j-1));
            cond2(j) = sum(wtstmpsort1(j+1:end));
        elseif j ==1
            cond2(j) = sum(wtstmpsort1(j+1:end));
        elseif j == length(kvectmp)
            cond1(j) = sum(wtstmpsort1(1:j-1));
        end
    end
    cond1 = cond1 <= wtstot/2;
    cond2 = cond2 <= wtstot/2;
    iswtdmed = and(cond1,cond2);
    if sum(iswtdmed) == 1
        wtdmedk(i) = kvectmpsort(iswtdmed);
    elseif sum(iswtdmed) == 0
        findcond1 = find(cond1);
        findcond2 = find(cond2);
        wtdmedk(i) = mean([kvectmpsort(findcond1(end)),kvectmpsort(findcond2(1))]);
    end

    medk(i) = median(kvectmp);
    sek(i) = std(kvectmp)/sqrt(numel(kvectmp));
    fplot(i) = median(fvectmp);
    num_bin(i) = numel(kvectmp);
end



end