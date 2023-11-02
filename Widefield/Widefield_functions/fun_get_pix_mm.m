% fun_get_pix_mm.m

function [pix_mm] = fun_get_pix_mm(recprams)

try
    mag = str2num(recprams.magnification)
catch
    paramsfile = dir('*Parameters.mat');
    paramstmp = load(paramsfile(1).name);
    mag = str2num(paramstmp.recprams.magnification)
end

if mag == 1
    pix_mm = 92.8
elseif mag == 1.25
    pix_mm = 115.0
elseif mag == 1.6
    pix_mm = 147.2
elseif mag == 2
    pix_mm = 183.2
elseif mag == 3.2
    pix_mm = 299.7
elseif mag == 4
    pix_mm = 377.3
elseif mag == 8
    pix_mm = 739.2
end

end