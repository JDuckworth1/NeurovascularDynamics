% NodeAngleFilter.m

function [angles] = NodeAngleFilter(nodetojoin,linkstojoin,vsl_graph_tmp,im_mask)
    
    masksize = size(im_mask);
    angles = zeros(max(size(nodetojoin)),1);
    for i=1:max(size(nodetojoin))
        ves1 = vsl_graph_tmp.link.cc_ind{linkstojoin(i,1),1}; %inds of vessels
        ves2 = vsl_graph_tmp.link.cc_ind{linkstojoin(i,2),1};
    
        [subr1,subc1] = ind2sub(masksize,ves1); %subs of vessels
        [subr2,subc2] = ind2sub(masksize,ves2);
        [noder,nodec] = ind2sub(masksize,vsl_graph_tmp.node.cc_ind{nodetojoin(i),1}); %subs of node
        mnoder = mean(noder); mnodec = mean(nodec);

        %Calc distance from node to each point in both vessels
        ndist1 = pdist2([mnodec,mnodec],[subc1,subr1]);
        if ndist1(end)<ndist1(1)
            ves1 = flip(ves1);
            [subr1,subc1] = ind2sub(masksize,ves1);
        end
        ndist2 = pdist2([mnodec,mnodec],[subc2,subr2]);
        if ndist2(end)<ndist2(1)
            ves2 = flip(ves2);
            [subr2,subc2] = ind2sub(masksize,ves2);
        end
        u1 = [mean(subc1)-mnodec,mean(subr1)-mnoder]; %vectors from node to each vessel link
        u2 = [mean(subc2)-mnodec,mean(subr2)-mnoder];

        angles(i) = acos(dot(u1,u2)/(norm(u1)*norm(u2))); %Angle between vessels, at joining node
    end

end
    




