%VesselDirection.m


function [ves_ang] = VesselDirection(lrs,is_near_midline_Q,rim)
    
    segs_link = lrs(1).segs_link;
    rim = rim.rim;
    ves_ang = zeros(length(segs_link),5);
    %Get angle of midline
    [midrow,midcol] = ind2sub(size(rim),find(is_near_midline_Q));

    midang = atan2(midrow(end)-midrow(1),midcol(end)-midcol(1));
    vang = zeros(length(segs_link),1);
    for i = 1:length(segs_link)
        pvec = zeros(segs_link(i),2);
%         indvec = zeros(segs_link(i),1);
        for j = 1:segs_link(i)
            pvec(j,1) = lrs(j).locs(i,1);
            pvec(j,2) = lrs(j).locs(i,2);
        end
%         indvec = sub2ind(size(im_mask),pvec(:,1),pvec(:,2));
        
        %Calc vector from first to last vessel pixel
        vvec = pvec(end,:) - pvec(1,:); % (rows,cols) (y,x)
        vvec = flip(vvec); % (cols,rows) (x,y)
        vang(i) = atan2(vvec(2),vvec(1)); %Radians
        vang_midline = vang(i)-midang;
        if vang_midline > pi
            vang_midline = vang_midline - 2*pi;
        elseif vang_midline < -pi
            vang_midline = vang_midline + 2*pi;
        end
        ves_ang(i,1) = vang_midline;
        ves_ang(i,2) = mean(pvec(:,1)); %Average row (y)
        ves_ang(i,3) = mean(pvec(:,2)); %Average col (x)

    end

    %Calculate location for each vessel in normalized coordinate system
    %Find midpoint
    midfit = polyfit(midcol,midrow,1); %Slope is midfit(1)
    midvals = polyval(midfit,1:size(rim,2)); %Fitted midline values
    midlineinds = sub2ind(size(rim),round(midvals),1:size(rim,2));
    midline_rim = intersect(find(rim),midlineinds); %Get part of midline that intersects the rim
    midpoint_ind = midline_rim(ceil(end/2)); %Midpoint of rim
    [midpointrow,midpointcol] = ind2sub(size(rim),midpoint_ind);

    %Now rotate everything around midpoint
    gam = -midang; %Make negative because y axis if flipped
    R = [cos(gam),-sin(gam);sin(gam),cos(gam)];
    %Get all x,y values to rotate
    centroids(:,1) = ves_ang(:,3); %x
    centroids(:,2) = ves_ang(:,2); %y 
    midptx = ones(length(centroids),1)*midpointcol;
    midpty = ones(length(centroids),1)*midpointrow;
    centroids_midpt = centroids - [midptx,midpty]; %(x,y)
    rotcentroids = R*centroids_midpt';
    rotcentroids = rotcentroids' + [midptx,midpty];

    [rimr,rimc] = ind2sub(size(rim),find(rim));
    midptx2 = ones(length(rimc),1)*midpointcol;
    midpty2 = ones(length(rimr),1)*midpointrow;
    rim_midpt = [rimc,rimr] - [midptx2,midpty2];
    rotrim = R*rim_midpt';
    rotrim = rotrim' + [midptx2,midpty2];

    %Find distance from rim edges to midpoint
    d1 = abs(max(rotrim(:,1))-midpointcol); %Dist from mid to right edge
    d2 = abs(min(rotrim(:,1))-midpointcol); %Dist from mid to left edge
    d3 = abs(max(rotrim(:,2))-midpointrow); %Dist from mid to top (top is actually left hemisphere here)
    d4 = abs(min(rotrim(:,2))-midpointrow); %Dist from mid to bottom

    %Normalize all centroid points and save
    normcentroids = rotcentroids - [midptx,midpty];
    ctrx = normcentroids(:,1);
    ctry = normcentroids(:,2);
    findright = ctrx>=0;
    findleft = ctrx<0;
    ctrx(findright) = ctrx(findright)/d1;
    ctrx(findleft) = ctrx(findleft)/d2;
    normcentroids(:,1) = ctrx;

    findup = ctry>=0;
    finddown = ctry<0;
    ctry(findup) = ctry(findup)/d3;
    ctry(finddown) = ctry(finddown)/d4;
    normcentroids(:,2) = ctry;

    ves_ang(:,4) = normcentroids(:,1);
    ves_ang(:,5) = normcentroids(:,2);
end