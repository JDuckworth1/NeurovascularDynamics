%WeightedAvg_SEM.m
%Uses an effective degrees of freedom

function [wtdavg,wtdSEM,wtdSD] = WeightedAvg_SEM(data,wts)

wtdavg = sum(data.*wts)/sum(wts);

n = numel(data);
wtdx2 = sum(wts.*data.^2)/sum(wts);
wtd2mom = wtdx2 - wtdavg^2;
df_factor = (sum(wts)^2)/((sum(wts)^2) - sum(wts.^2));
df_a = (sum(wts)^2);
df_b = sum(wts.^2);
wtdVar = wtd2mom*df_factor;
wtdSD = sqrt(wtdVar);
wtdSEM = sqrt(wtd2mom*(df_b/(df_a-df_b)));

end