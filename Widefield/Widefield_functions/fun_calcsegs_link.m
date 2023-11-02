% fun_calcsegs_link.m

function [segs_link] = fun_calcsegs_link(link_results_struct)

        test = zeros(length(link_results_struct),length(link_results_struct(1).vsl_graph.link.cc_ind));
        for k=1:length(link_results_struct)
            for j=1:length(link_results_struct(1).vsl_graph.link.cc_ind)
                test(k,j) = link_results_struct(k).wave(1,j);  %First wave value from all segments, all links (check for NaNs)
            end
        end
        segs_link = zeros(length(link_results_struct(1).vsl_graph.link.cc_ind),1);
        for j=1:length(link_results_struct(1).vsl_graph.link.cc_ind)
            segments_vec = test(:,j);
            non_nanstmp = ~isnan(segments_vec);
            segs_link(j) = sum(non_nanstmp); %Number of segments in each link
        end

end