%WeightedAvg_SEM.m

function [wtdavg,wtdSEM,wtdSD] = WeightedAvg_SEM(data,wts)

wtdavg = sum(data.*wts)/sum(wts);
n = numel(data);
wtdx2 = sum(wts.*data.^2)/sum(wts);
wtd2mom = wtdx2 - wtdavg^2;
wtdVar = wtd2mom*n/(n-1);
wtdSD = sqrt(wtdVar);
wtdSEM = wtdSD/sqrt(n);

end