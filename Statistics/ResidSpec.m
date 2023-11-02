%ResidSpec.m
%Uses ftestc, rmlinesc from chronux
% http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/continuous/ftestc.html
% [Fval,A,f,sig,sd] = ftestc(data,params,p,plt)

function [data_nolines,f,Sresid,Stot,Amps,fmax] = ResidSpec(data,params)

p = 0.05;
p = p/length(data);
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear err trialave
data=change_row_to_column(data);
[N,C]=size(data);
tapers=dpsschk(tapers,N,Fs); % calculate the tapers
[N,K]=size(tapers);
nfft=max(2^(nextpow2(N)+pad),N);% number of points in fft
[f,findx]=getfgrid(Fs,nfft,fpass);% frequency grid to be returned

Kodd=1:2:K;
Keven=2:2:K;
J=mtfftc(data,tapers,nfft,Fs);% tapered fft of data - f x K x C
Jp=J(findx,Kodd,:); % drop the even ffts and restrict fft to specified frequency grid - f x K x C
tapers=tapers(:,:,ones(1,C)); % add channel indices to the tapers - t x K x C
H0 = squeeze(sum(tapers(:,Kodd,:),1)); % calculate sum of tapers for even prolates - K x C
if C==1;H0=H0';end

Nf=length(findx);% number of frequencies
H0 = H0(:,:,ones(1,Nf)); % add frequency indices to H0 - K x C x f
H0=permute(H0,[3 1 2]); % permute H0 to get dimensions to match those of Jp - f x K x C
H0sq=sum(H0.*H0,2);% sum of squares of H0^2 across taper indices - f x C
JpH0=sum(Jp.*squeeze(H0),2);% sum of the product of Jp and H0 across taper indices - f x C
A=squeeze(JpH0./H0sq); % amplitudes for all frequencies and channels
Kp=size(Jp,2); % number of even prolates
Ap=A(:,:,ones(1,Kp)); % add the taper index to C
Ap=permute(Ap,[1 3 2]); % permute indices to match those of H0
Jhat=Ap.*H0; % fitted value for the fft

num=(K-1).*(abs(A).^2).*squeeze(H0sq);%numerator for F-statistic
den=squeeze(sum(abs(Jp-Jhat).^2,2)+sum(abs(J(findx,Keven,:)).^2,2));% denominator for F-statistic
Fval=num./den; % F-statisitic
sig=finv(1-p,2,2*K-2); % F-distribution based 1-p% point
A=A*Fs;

%From fitlinesc.m
fmax=findpeaks(Fval,sig);
freqs=cell(1,C);
Amps=cell(1,C);
datafit=data;
for ch=1:C
    fsig=f(fmax(ch).loc);
    freqs{ch}=fsig;
    Amps{ch}=A(fmax(ch).loc,ch);
    Nf=length(fsig);
    datafit(:,ch)=exp(i*2*pi*(0:N-1)'*fsig/Fs)*A(fmax(ch).loc,ch)+exp(-i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(fmax(ch).loc,ch)); %These are the sinewaves
end

data_nolines = data - datafit;
[Stot,f]=mtspectrumc(data,params);
[Sresid,f]=mtspectrumc(data_nolines,params);

end




