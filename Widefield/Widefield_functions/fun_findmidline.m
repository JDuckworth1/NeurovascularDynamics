% fun_findmidline.m

function [is_near_midline_Q] = fun_findmidline(animal,im_mask,rim,toplot)

if ~strcmp(animal,'JD221024F2') && ~strcmp(animal,'JD221024F3') && ~strcmp(animal,'TB200217F2')
    try
        [is_near_midline_Q,mix_gaussian_model] = Thomas_generate_hemisphere_mask(toplot.mask,uint8(im_mask),rim);
        assert(abs(diff(mix_gaussian_model.mu(:, 2))) > 350)
    catch
        try
            [is_near_midline_Q,mix_gaussian_model] = Thomas_generate_hemisphere_mask(toplot.mask,uint8(im_mask),rim);
            assert(abs(diff(mix_gaussian_model.mu(:, 2))) > 350)
        catch
            [is_near_midline_Q,mix_gaussian_model] = Thomas_generate_hemisphere_mask(toplot.mask,uint8(im_mask),rim);
            assert(abs(diff(mix_gaussian_model.mu(:, 2))) > 350)
        end
    end

end


end