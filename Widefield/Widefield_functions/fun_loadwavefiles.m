% fun_loadwavefiles.m
%Loads the current trial's wave (df/f for vessels and neurons, as well as
%experimental parameters)

function [wave,neurowave,toplot,neurotoplot,recprams] = fun_loadwavefiles(data_folder,filename,stim)

if stim == 1
    try
        wave = h5read(append(data_folder,'\',filename,'svd_puff_wave.h5'),'/wave');
    catch
        try
            wave = h5read(append(data_folder,'\',filename,'puff_wave.h5'),'/wave');
        catch
            wave = h5read(append(data_folder,'\',filename,'_puff_wave.h5'),'/wave');
        end
    end
    try
        neurowave = h5read(append(data_folder,'\',filename,'svd_puff_neurons_wave.h5'),'/wave');
    catch
        try
            neurowave = h5read(append(data_folder,'\',filename,'puff_neurons_wave.h5'),'/wave');
        catch
            neurowave = h5read(append(data_folder,'\',filename,'_puff_neurons_wave.h5'),'/wave');
        end
    end
    try
        toplottmp = load(append(data_folder,'\',filename,'svd_puff_toplot.mat'));
        toplot = toplottmp.toplot;
        recprams = toplottmp.recprams;
    catch
        try
            toplottmp = load(append(data_folder,'\',filename,'svd_puff.mat'));
            toplot = toplottmp.toplot;
            recprams = toplottmp.recprams;
        catch
            toplottmp = load(append(data_folder,'\',filename,'_puff.mat'));
            toplot = toplottmp.toplot;
            recprams = toplottmp.recprams;
        end
    end
    try
        neurotoplot = load(append(data_folder,'\',filename,'svd_puff_neurons_toplot.mat'));
        neurotoplot = neurotoplot.toplot;
    catch
        try
            tmpneurotoplot = load(append(data_folder,'\',filename,'svd_puff_neurons.mat'));
            neurotoplot = tmpneurotoplot.toplot;
        catch
            tmpneurotoplot = load(append(data_folder,'\',filename,'_puff_neurons.mat'));
            neurotoplot = tmpneurotoplot.toplot;
        end
    end
else
    try
        wave = h5read(append(data_folder,'\',filename,'svd_beforepuff_wave.h5'),'/wave');
    catch
        wave = h5read(append(data_folder,'\',filename,'_beforepuff_wave.h5'),'/wave');
    end
    try
        neurowave = h5read(append(data_folder,'\',filename,'svd_beforepuff_neurons_wave.h5'),'/wave');
    catch
        neurowave = h5read(append(data_folder,'\',filename,'_beforepuff_neurons_wave.h5'),'/wave');
    end
    try
        toplottmp = load(append(data_folder,'\',filename,'svd_beforepuff_toplot.mat'));
        toplot = toplottmp.toplot;
        recprams = toplottmp.recprams;
    catch
        try
            toplottmp = load(append(data_folder,'\',filename,'svd_beforepuff.mat'));
            toplot = toplottmp.toplot;
            recprams = toplottmp.recprams;
        catch
            toplottmp = load(append(data_folder,'\',filename,'_beforepuff.mat'));
            toplot = toplottmp.toplot;
            recprams = toplottmp.recprams;
        end
    end
    try
        neurotoplot = load(append(data_folder,'\',filename,'svd_beforepuff_neurons_toplot.mat'));
        neurotoplot = neurotoplot.toplot;
    catch
        try
            tmpneurotoplot = load(append(data_folder,'\',filename,'svd_beforepuff_neurons.mat'));
            neurotoplot = tmpneurotoplot.toplot;
        catch
            tmpneurotoplot = load(append(data_folder,'\',filename,'_beforepuff_neurons.mat'));
            neurotoplot = tmpneurotoplot.toplot;
        end
    end
end

end




