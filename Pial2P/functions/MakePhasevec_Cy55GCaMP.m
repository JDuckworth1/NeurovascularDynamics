

function [phasevec] = MakePhasevec_Cy55GCaMP(distvec,phivec,f_peak)

fit2 = polyfit(distvec,phivec,1);
yhat = polyval(fit2,distvec);
R = corrcoef(distvec,phivec); %Alternative is to use R2 formula, is the same as corrcoef for linear fits
SSE = sum((phivec-yhat).^2);
xbar = sum(distvec)/numel(distvec);
Sxx = sum((distvec-xbar).^2);
sigma = sqrt(numel(distvec)*SSE/((numel(distvec)-2)*Sxx));

rad_mm = fit2(1);
phasevec(1) = rad_mm; %rad_pix*pix_mm = k rad/mm
phasevec(2) = R(1,2);
phasevec(3) = f_peak;
phasevec(4) = sigma; %STD of slope (phase gradient)
phasevec(5) = distvec(end); %Length of vessel segment
phasevec(6) = sqrt(SSE/(numel(distvec)-2)); %Alternative measure of goodness of fit to R2
phasevec(7) = 2*pi/rad_mm*f_peak; 
phasevec(8) = fit2(2);

end
