%weightedMed_SinglePoint.m

function [wtdmed] = weightedMed_SinglePoint(kvec,wts)

wtstot = sum(wts);
[kvecsort,sortinds] = sort(kvec); %Sort k in increasing order
wtstmpsort1 = wts(sortinds); %Sort weights by k order

cond1 = zeros(length(kvec),1);
cond2 = zeros(length(kvec),1);
for j = 1:length(kvec)
    cond1(j) = sum(wtstmpsort1(1:j-1));
    cond2(j) = sum(wtstmpsort1(j+1:end));
end
cond1 = cond1 <= wtstot/2;
cond2 = cond2 <= wtstot/2;
iswtdmed = and(cond1,cond2);
if sum(iswtdmed) == 1
    wtdmed = kvecsort(iswtdmed);
elseif sum(iswtdmed) == 0
    findcond1 = find(cond1);
    findcond2 = find(cond2);
    wtdmed = mean([kvecsort(findcond1(end)),kvecsort(findcond2(1))]);
end


end