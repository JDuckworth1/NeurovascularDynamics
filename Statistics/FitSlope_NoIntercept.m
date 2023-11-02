% FitSlope_NoIntercept.m

% Use least-squares formulation to estimate slope b from data for best fit
% line y = b*x. Also calcualate the uncertainty in b, from uncertainties in
% each point y (weights). 

function [b,SEb_origin,R2Origin] = FitSlope_NoIntercept(xpts,ypts)

    bnum = sum(xpts.*ypts);
    bdenom = sum(xpts.^2);
    b = bnum/bdenom; %Least-squares slope
    n = numel(xpts);

    %Calculate slope uncertainty
    %This is using the formula directly from linear regression equations
    %(y = mx+b)

%     meany = mean(ypts);
%     meanx = mean(xpts);
% 
%     totSS = sum((ypts-meany).^2);
%     regSS = (sum((xpts-meanx).*(ypts-meany))^2)/sum((xpts-meanx).^2);
%     residMS = (totSS - regSS)/(n-1);
%     Varb = residMS/sum(xpts.^2);
%     SEb = sqrt(Varb)

    clearvars -except xpts ypts meanx meany n b SEb
    %This is from regression through the origin equations, using 0 as
    %reference point instead mean(y) 
    %This is also how matlab calculates SEb for CIs when forcing int=0
    totSS_o = sum(ypts.^2);
    regSS = (sum(xpts.*ypts)^2)/sum(xpts.^2);
    residMS = (totSS_o - regSS)/(n-1);
    Varb = residMS/sum(xpts.^2);
    SEb_origin = sqrt(Varb);

    %Calculate R2 using alternative definition
    %https://web.ist.utl.pt/~ist11038/compute/errtheory/,regression/regrthroughorigin.pdf
    yhat = xpts*b;
    SSR = sum(yhat.^2);
    SST = sum(ypts.^2);
    R2Origin = SSR/SST;

    top = sum((yhat-mean(ypts)).^2);
    bot = sum((ypts-mean(ypts)).^2);
    normalR2 = top/bot;

end


