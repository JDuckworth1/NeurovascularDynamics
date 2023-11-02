% fun_get_PAMag_depth.m

function [PA,depth1,depth2,zoom,pix_um] = fun_get_PAMag_depth(filename)

    if contains(filename,'6.0x')
        namestr = extractBefore(filename,'_6.0x');
    elseif contains(filename,'7.6x')
        namestr = extractBefore(filename,'_7.6');
    elseif contains(filename,'_6x')
        namestr = extractBefore(filename,'_6x');
    elseif contains(filename,'_6.0x')
        namestr = extractBefore(filename,'_6.0x');
    elseif contains(filename,'6.4x')
        namestr = extractBefore(filename,'_6.4x');
    elseif contains(filename,'3.0x')
        namestr = extractBefore(filename,'_3.0x');
    elseif contains(filename,'5.0x')
        namestr = extractBefore(filename,'_5.0x');
    elseif contains(filename,'5.5x')
        namestr = extractBefore(filename,'_5.5x');
    elseif contains(filename,'5.6x')
        namestr = extractBefore(filename,'_5.6x');
    elseif contains(filename,'5.8x')
        namestr = extractBefore(filename,'_5.8x');
    elseif contains(filename,'5.7x')
        namestr = extractBefore(filename,'_5.7x');
    elseif contains(filename,'7.0x')
        namestr = extractBefore(filename,'_7.0x');
    elseif contains(filename,'7x')
        namestr = extractBefore(filename,'_7x');
    elseif contains(filename,'6.5x')
        namestr = extractBefore(filename,'_6.5x');
    elseif contains(filename,'4.5x')
        namestr = extractBefore(filename,'_4.5x');
    elseif contains(filename,'4.8x')
        namestr = extractBefore(filename,'_4.8x');
    else
        disp('Magnification not found')
        return
    end

    PA = extractBefore(namestr,'_');
    PA = extractAfter(PA,2)
    
    depth_str = extractAfter(namestr,[PA,'_']);
    depth1 = extractBefore(depth_str,'_')
    depth2 = extractAfter(depth_str,'_')

    try
        zoom = extractBetween(filename,[depth2,'_'],'x_00001.tif');
    catch
        try
            zoom = extractBetween(filename,[depth2,'_'],'x_stim_00001.tif');
        catch
            zoom = extractBetween(filename,[depth2,'_'],'x_stim_rest_00001.tif');
        end
    end
    zoom = str2double(zoom{1})

    pix_um = 300/100;
    pix_um = pix_um / 7.6 * zoom

end




