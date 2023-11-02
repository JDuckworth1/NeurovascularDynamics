% fun_get_NrOCh.m

function [NrOCh] = fun_get_NrOCh(str,fname,str1,str2)

if contains(fname,str1);
    NrOCh = '2G' %options 1 = one channel, 2G = green of dual, 2R = red, 2 for Ratiometricmeasurements. of dual this is inversed when using SMART to control the LED exposure
    if contains(fname,'ASAP') && str =='v';
        NrOCh = '2R' %options 1 = one channel, 2G = green of dual, 2R = red, 2 for Ratiometricmeasurements. of dual this is inversed when using SMART to control the LED exposure
    end
elseif contains(fname,str2);
    NrOCh='1'
else
    NrOCh = input('How many channels where imageged 2G, 1, 2')
end



end




