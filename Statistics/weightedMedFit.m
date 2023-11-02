

function [wtdmedf,wtdmedk] = weightedMedFit(kvec,fvec,wts,numbins)


% Fit f vs k 10% f bins, weighted mean k in each bin, weighted median k
[f,d] = ecdf(fvec);

binvec = zeros(numbins,1);
for i=1:numbins
    [~,isclosestind] = min(abs(f-i/numbins));
    binvec(i) = d(isclosestind);
end
binvec = [0;binvec];
fbins = binvec';
meank = zeros(size(fbins,2)-1,1);
wtdmedk = zeros(size(fbins,2)-1,1);
wtdmedf = zeros(size(fbins,2)-1,1);
sek = zeros(size(fbins,2)-1,1);
fplot = zeros(size(fbins,2)-1,1);
num_bin = zeros(size(fbins,2)-1,1);

for i=2:size(fbins,2)
    %Get f, k, weights in current bin
    inbin1 = fvec > fbins(i-1);
    if i==size(fbins,2)
        inbin2 = fvec <= max(fvec);
    else
        inbin2 = fvec <= fbins(i);
    end
    inbin = and(inbin1,inbin2);
    fvectmp = fvec(inbin);
    kvectmp = kvec(inbin);
    wtstmp = wts(inbin);

    meannumk = dot(wtstmp,kvectmp);
    meandenom = sum(wtstmp);
    meank(i-1) = meannumk/meandenom; %Weighted mean k within bin

    wtstot = sum(wtstmp);
    [kvectmpsort,sortinds] = sort(kvectmp); %Sort k in increasing order
    wtstmpsort1 = wtstmp(sortinds); %Sort weights by k order

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
        wtdmedk(i-1) = kvectmpsort(iswtdmed);
    elseif sum(iswtdmed) == 0
        findcond1 = find(cond1);
        findcond2 = find(cond2);
        wtdmedk(i-1) = mean([kvectmpsort(findcond1(end)),kvectmpsort(findcond2(1))]);
    end

    %Calc weighted f for plotting
    [fvectmpsort2,sortinds2] = sort(fvectmp); %Sort f by increasing order
    wtstmpsort2 = wtstmp(sortinds2); %Weights sorted by f order

    cond1 = zeros(length(fvectmp),1);
    cond2 = zeros(length(fvectmp),1);
    for j = 1:length(fvectmp)
        if j ~= 1 && j ~= length(fvectmp)
            cond1(j) = sum(wtstmpsort2(1:j-1));
            cond2(j) = sum(wtstmpsort2(j+1:end));
        elseif j ==1
            cond2(j) = sum(wtstmpsort2(j+1:end));
        elseif j == length(fvectmp)
            cond1(j) = sum(wtstmpsort2(1:j-1));
        end
    end
    cond1 = cond1 <= wtstot/2;
    cond2 = cond2 <= wtstot/2;
    iswtdmed = and(cond1,cond2);
    if sum(iswtdmed) == 1
        wtdmedf(i-1) = fvectmpsort2(iswtdmed);
    elseif sum(iswtdmed) == 0
        findcond1 = find(cond1);
        findcond2 = find(cond2);
        wtdmedf(i-1) = mean([fvectmpsort2(findcond1(end)),fvectmpsort2(findcond2(1))]);
    end

%     medk(i-1) = median(kvectmp);
%     sek(i-1) = std(kvectmp)/sqrt(numel(kvectmp));
%     fplot(i-1) = median(fvectmp);
%     num_bin(i-1) = numel(kvectmp);
end

end