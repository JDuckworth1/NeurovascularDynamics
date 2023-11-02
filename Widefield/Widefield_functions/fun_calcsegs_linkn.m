% fun_calcsegs_link.m

function [segs_linkn] = fun_calcsegs_link(link_results_struct_neu)

        testn = zeros(length(link_results_struct_neu),length(link_results_struct_neu(1).n_vsl_graph.link.cc_ind));
        for k=1:length(link_results_struct_neu)
            for j=1:length(link_results_struct_neu(1).n_vsl_graph.link.cc_ind)
                testn(k,j) = link_results_struct_neu(k).wave(1,j);  %First wave value from all segments, all links (check for NaNs)
            end
        end
        segs_linkn = zeros(length(link_results_struct_neu(1).n_vsl_graph.link.cc_ind),1);
        for j=1:length(link_results_struct_neu(1).n_vsl_graph.link.cc_ind)
            segments_vecn = testn(:,j);
            non_nanstmpn = ~isnan(segments_vecn);
            segs_linkn(j) = sum(non_nanstmpn); %Number of segments in each link
        end

end