% widefield_ExtractVesNeuPhaseGrad.m

% % % % % 

% Results for each trial are saved in "link_results_struct" with each row
% representing consecuative vessel segments within each vessel link.  Entries
% include:
% wave: df/f (frames x links). A link with 10 segments will have 10 "wave"
%   rows, with the 11th wave row being all NaN values.
% locs: location of each segment within the imaging frame (row,column).
% linknum: pixels along the vessel skeleton required to reach each segment
%   (from the initial segment).
% skel_lab: skeleton segment number, as used in toplot.skel_label.
% Other results are saved in link_results_struct row 1, and are later combined
% across trials to produce ensemble results.

% "Segment" refers to a group of pixels within each vessel / neuron
% envelope. "Link" refers to the part of the vessel network that spans two
% nodes. The vessel network is divided into links, which are then filtered
% by length for further analysis.

% Contact: Jacob Duckworth, jaduckwo@ucsd.edu
% % % % %


clc;clear;close all;
animal = 'TB013120M2'
files = dir(['*',animal,'*\*meters.mat'])
stim = 0; %Analyze stim period? (0 for before-stim period).
savefolder = 'AllSegments\'; %Choose where to save

%% Iterates over all trials for a given animal
for file = 1:length(files)
close all
clearvars -except files file animal stim

data_folder = files(file).folder;
filename = extractBefore(files(file).name,'_Parameters.mat');

[folder_name] = fun_get_foldername(filename);
if stim == 1
    cd([savefolder,'Stim\',animal])
else
    cd([savefolder,animal])
end

mkdir(folder_name)
if stim == 1
    savefile = append(savefolder,'Stim\',animal,'\',folder_name);
elseif stim == 0
    savefile = append(savefolder,animal,'\',folder_name);
end

%Load dF/F files and plotting files
cd(data_folder)
[wave,neurowave,toplot,neurotoplot,recprams] = fun_loadwavefiles(data_folder,filename,stim);

try
    im_mask = imread('Mask.tif');
catch
    im_mask = toplot.mask;
end
rim = load('rim.mat');

%Magnification recorded with trial info
[pix_mm] = fun_get_pix_mm(recprams);
%Create vessel graph for current trial
[toplot,vsl_graph,im_mask] = fun_createvesgraph(im_mask,toplot);
%Get window midline with mix gaussian model
[is_near_midline_Q] = fun_findmidline(animal,im_mask,rim,toplot);

modifygraph = 1;
%% Modify the graph to exclude small branching segments and join the remaining ones
if modifygraph == 1
    map = zeros(toplot.mask_size);
    map(toplot.mask_ind) = 1;
    map(vsl_graph.link.pos_ind) = -1;
    map(vsl_graph.node.pos_ind) = 2;
    figure; subplot(1,2,1); imagesc(map); daspect([1,1,1]); ylim([550 1150]); xlim([50 550]); axis off;
    title('Original Graph')
    
%% --------------------- Step 1: delete small links that are only connected to one node ---------------------
    minvox = 0.4*pix_mm %Size threshold for small links (pixels)
    vox_per_cc_tmp = vsl_graph.link.num_voxel_per_cc; %Pixels per each link
    issmall = vox_per_cc_tmp < minvox;
    isisolated = vsl_graph.link.num_node < 2;
    todel = and(issmall,isisolated);

    vsl_graph_tmp = vsl_graph; %Make copy of vsl_graph
    links = vsl_graph_tmp.link.cc_ind; %Linear indices of all links
    links(todel) = []; %Delete short links that only have one node
    vsl_graph_tmp.link.cc_ind = []; %Replace cc_inds after deleting 
    vsl_graph_tmp.link.cc_ind = links;
    vsl_graph_tmp.link.pos_ind = []; %Replace list of all link indices
    vsl_graph_tmp.link.pos_ind = cell2mat(vsl_graph_tmp.link.cc_ind(:)); %Concatenate remaining link indices
    cnltmp = vsl_graph_tmp.link.connected_node_label;
    cnltmp(todel,:) = []; %Delete node labels for deleted links
    vsl_graph_tmp.link.connected_node_label = cnltmp;
    vsl_graph_tmp.link.num_voxel_per_cc(todel) = []; %Delete link lengths for deleted links
    numnodetmp = vsl_graph_tmp.link.num_node;
    numnodetmp(todel) = []; %Replace number of nodes for each remaining link
    vsl_graph_tmp.link.num_node = numnodetmp;

    vsl_graph_tmp.link.num_cc = sum(~todel); %Remaining links in the graph
    vsl_graph_tmp.link.num_voxel = size(vsl_graph_tmp.link.pos_ind,1);
    vsl_graph_tmp.link.label = repelem(1:vsl_graph_tmp.link.num_cc, vsl_graph_tmp.link.num_voxel_per_cc)';
    vsl_graph_tmp.link.map_ind_2_label = sparse(vsl_graph_tmp.link.pos_ind, ...
        ones(vsl_graph_tmp.link.num_voxel,1), ...
        vsl_graph_tmp.link.label, ...
        vsl_graph_tmp.num.block_voxel,1);
    %Now all "link" fields are updated as necessary

    map = zeros(size(toplot.mask));
    map(toplot.mask_ind) = 1;
    map(vsl_graph_tmp.link.pos_ind) = -1;
    map(vsl_graph_tmp.node.pos_ind) = 2;
    subplot(1,2,2); imagesc(map); daspect([1,1,1]); ylim([550 1150]); xlim([50 550]); axis off;
    str = sprintf('Min link size = %.2f mm',minvox/pix_mm);
    title(str) %Plot graph after eliminating small segments
    
%% --------------------- Step 2: Join together links across nodes ---------------------
    counter = 1;
    nodetojoin = 0;
    for i=1:length(links) %Iterate over number of remaining links
        node(1) = vsl_graph_tmp.link.connected_node_label(i,1); %find which nodes are connected to each segment
        node(2) = vsl_graph_tmp.link.connected_node_label(i,2); %find which nodes are connected to each segment
        findnode1 = find(vsl_graph_tmp.link.connected_node_label==node(1)); %Find where else this node is connected
        findnode2 = find(vsl_graph_tmp.link.connected_node_label==node(2)); %Find where else this node is connected
        
        if length(findnode1) == 2 && ~ismember(node(1),nodetojoin) && all(node) %Create list of nodes to join connected links around
            nodetojoin(counter) = node(1);
            counter = counter +1;
        end
        if length(findnode2) == 2 && ~ismember(node(2),nodetojoin) && all(node)
            nodetojoin(counter) = node(2);
            counter = counter + 1;
        end
    end

    %Identify links to be joined at each node in nodetojoin
    linkstojoin = zeros(size(nodetojoin,2),2);
    for i=1:length(nodetojoin) %Iterate over number of nodes to be joined
        %Make linkstojoin matrix: nodes x 2 links that are being joined
        [linkstojoin(i,:),~] = ind2sub(size(vsl_graph_tmp.link.connected_node_label),find(vsl_graph_tmp.link.connected_node_label==nodetojoin(i)));
        %The row in conn_node_label is which link. There are 2 for each node in node2join-> linkstojoin is nodes x 2
    end

    %Use link angle threshold function
    anglethresh = pi/2;
    [angles] = NodeAngleFilter(nodetojoin,linkstojoin,vsl_graph_tmp,im_mask); %Angle across node between each links to be joined
    nodetojoin(angles<anglethresh) = [];
    linkstojoin(angles<anglethresh,:) = []; %Delete nodes with angles smaller than the threshold from nodeto join and linkstojoin. (Don't want to join nodes with small angles, like MCA big branch)

    %Gather all links that are connected one-after the other
    vsl_graph_tmp.link.seglinks = zeros(length(linkstojoin),100);
    vsl_graph_tmp.link.segrows = zeros(length(linkstojoin),100);
    for i = 1:size(linkstojoin,1)
        repeat = 1;
        segrows(1) = i;
        seglinks(1) = linkstojoin(i,1);
        counter = 2;
        whilecounter = 1;
        circletest = 0;
        while repeat == 1
            if whilecounter == 1 %For the first link
                linkstmp = linkstojoin(i,:);
                entrytmp = linkstmp(1);
                entryrow = i;
            else
                entrytmp = newentry;
                entryrow = newrow;
            end
            test = find(linkstojoin == entrytmp);
            if length(test)==2 && circletest == 0
                [rows,~] = ind2sub(size(linkstojoin),test'); %Get rows where the entry is found (row of link2join is node)
                newrow = rows(rows~=entryrow); %Save this | Get the new row (next one to connect)
                newentries = linkstojoin(newrow,:);
                newentry = newentries(newentries~=entrytmp); %New entry (link) is in the new row, not the old link. Now repeat
                
                segrows(counter) = newrow;
                seglinks(counter) = newentry;
                counter = counter +1;
                whilecounter = whilecounter +1;
                if newentry == linkstmp(1)
                    circletest = 1;
                end
            else
                repeat = 0;
            end
        end
        %Repeat for second link in linkstojoin vector:
        repeat = 1;
        segrows2(1) = i;
        seglinks2(1) = linkstojoin(i,2);
        counter = 2;
        whilecounter = 1;
        circletest = 0;
        while repeat == 1
            if whilecounter == 1 %For the first link
                linkstmp = linkstojoin(i,:);
                entrytmp = linkstmp(2); %(Second column)
                entryrow = i;
            else
                entrytmp = newentry;
                entryrow = newrow;
            end
            test = find(linkstojoin == entrytmp);
            if length(test)==2 && circletest == 0 %See above was 1
                [rows,~] = ind2sub(size(linkstojoin),test'); %Get rows where the entry is found
                newrow = rows(rows~=entryrow); %Save this | this will fail when length(test) = 0
                newentries = linkstojoin(newrow,:);
                newentry = newentries(newentries~=entrytmp); %Now repeat
                
                segrows2(counter) = newrow;
                seglinks2(counter) = newentry;
                counter = counter +1;
                whilecounter = whilecounter +1;
                if newentry == linkstmp(1)
                    circletest = 1;
                end
            else
                repeat = 0;
            end
        end
        
        comb_seglinks = unique([seglinks,seglinks2],'stable'); %Combine found links from original link 1 and 2 (get rid of repeated links, don't artificially order them)
        comb_segrows = unique([segrows,segrows2],'stable');
        vsl_graph_tmp.link.seglinks(i,1:length(comb_seglinks)) = comb_seglinks;
        vsl_graph_tmp.link.segrows(i,1:length(comb_segrows)) = comb_segrows; %Remember rows are in rows of linkstojoin, not original connected_node_label
        
        clear segrows2 segrows seglinks2 seglinks
    end
    tmp_seglinks = vsl_graph_tmp.link.seglinks;
    tmp_segrows = vsl_graph_tmp.link.segrows;
    
    %Sort these Results and delete repeated connected links
    tmpmat = vsl_graph_tmp.link.seglinks;
    tmpmat2 = sort(tmpmat,2); %Need to do this do delete repeated rows, sorts each row in ascending order, puts all the zeros on the left
    tmpmat3 = unique(tmpmat2,'rows'); %Gives the unique rows, sorts them in ascending order
    tmpmatr = vsl_graph_tmp.link.segrows;
    tmpmatr2 = sort(tmpmatr,2); %Need to do this do delete repeated rows (sort each ROW in ascending order)
    tmpmatr3 = unique(tmpmatr2,'rows'); %Treats each row as entity, returns unique rows in sorted order, by #zeros and size of first entry in row
    vsl_graph_tmp.link.seglinks = tmpmat3;
    vsl_graph_tmp.link.segrows = tmpmatr3;
    
    %Connect these links and update vsl_graph_tmp accordingly
    for i=1:size(vsl_graph_tmp.link.seglinks,1) %Iterate over number of links to connect
        linkstmp = nonzeros(vsl_graph_tmp.link.seglinks(i,:));
        linktoupdate = linkstmp(1); %Update first link in group
        findlink = find(tmp_seglinks==linktoupdate);
        [tmp_seglinks_row,~] = ind2sub(size(tmp_seglinks),findlink(1)); %Get row of tmp_seglinks that linktoupdate first appears in (that row in tmp_segrows will give the corresponding nodes)
        nodes = nonzeros(tmp_segrows(tmp_seglinks_row,:)); %Want row of tmp_seglinks that linktoupdate appears in

        for j = 1:length(linkstmp)-1 
            nodeinds = vsl_graph_tmp.node.cc_ind{nodetojoin(nodes(j)),1}; %Get the 2-4 node inds corresponding to relevant link. segrows are which node in nodetojoin
            nextlink = vsl_graph_tmp.link.cc_ind{linkstmp(j+1),1};
            vsl_graph_tmp.link.cc_ind{linktoupdate,1} = [vsl_graph_tmp.link.cc_ind{linktoupdate,1};nextlink]; %Here we don't add nodes (will delete/add later)
        end %Need to delete second entry of linkstmp/seglinks because it's concatenated into the first entries indices.
        vsl_graph_tmp.link.cc_ind{linktoupdate,1} = unique(vsl_graph_tmp.link.cc_ind{linktoupdate,1},"stable"); %Get rid of repeated entries, should be ordered directionally already so don't change that
    end
    
    %% Step 2.2

    for seglinkiter = 1:size(vsl_graph_tmp.link.seglinks,1)
        seglinktmp = nonzeros(vsl_graph_tmp.link.seglinks(seglinkiter,:));
        startlink = seglinktmp(1);
        %Calc avg. loc of all links then order them accordingly.
        for i=1:length(seglinktmp)
            seglinklength(i) = length(vsl_graph_tmp.link.cc_ind{seglinktmp(i),1});
        end
        seglinkinds = NaN(length(seglinktmp),max(seglinklength));
        for i=1:length(seglinktmp) %get indices for each segment
            seglinkinds(i,1:seglinklength(i)) = vsl_graph_tmp.link.cc_ind{seglinktmp(i),1};
        end
        %Delete added entries from the first row
        seglinkindstmp1 = seglinkinds(1,:); %First row (has had entries concatenated onto it)
        seglinkindstmp2 = seglinkinds(2:end,:); %Rest of the rows (just raw segments)
        seglinkindstmp2 = seglinkindstmp2(~isnan(seglinkindstmp2));
        seglinkindstmp1 = seglinkindstmp1(~ismember(seglinkindstmp1,seglinkindstmp2)); %Delete added pixels
        seglinkinds(1,1:length(seglinkindstmp1)) = seglinkindstmp1; %Seglinkinds is now just each individual link
        seglinkinds(1,(length(seglinkindstmp1)+1):end) = NaN;
        %Add nodes to links
        findlink = find(tmp_seglinks==startlink);
        [tmp_seglinks_row,~] = ind2sub(size(tmp_seglinks),findlink(1));
        nodes = nonzeros(tmp_segrows(tmp_seglinks_row,:));
        
        nodeinds1 = nodetojoin(nodes);
        nodeinds = zeros(length(nodeinds1),10); %Just add extra zeros (no nodes are 10 pixels wide)
        for i=1:length(nodeinds1)
            nodeinds(i,1:length(vsl_graph_tmp.node.cc_ind{nodeinds1(i),1})) = vsl_graph_tmp.node.cc_ind{nodeinds1(i),1};
        end %Now delete these inds from seglinkinds and re-add them to appropriate links
        
        seglinkindstmp1 = seglinkindstmp1(~ismember(seglinkindstmp1,nodeinds)); %Delete nodes from first row of seglinkinds
        %Calc avg position of each node: Replaced prev. version 7/10
        for i=1:size(nodeinds,1)
           tmp_nodeinds = nodeinds(i,:);
           [row,col] = ind2sub(size(im_mask),nonzeros(tmp_nodeinds)); 
           nodepos(i,1) = mean(row);
           nodepos(i,2) = mean(col);
        end
        %Calc position of each segment endpoint
        for i=1:size(seglinkinds,1)
           tmp_linkinds = seglinkinds(i,:);
           tmp_linkinds = tmp_linkinds(~isnan(tmp_linkinds));
           [link_startpos(i,1),link_startpos(i,2)] = ind2sub(size(im_mask),tmp_linkinds(1));
           [link_endpos(i,1),link_endpos(i,2)] = ind2sub(size(im_mask),tmp_linkinds(end)); 
        end
        %Find which link pixel each node is closest to
        %For each node, need to calculate (2*links) distances
        nodedist = zeros(size(nodeinds,1),2*size(seglinkinds,1));
        linkiter = 1:2:2*size(seglinkinds,1);
        for node = 1:size(nodeinds,1)
            counter = 1;
            for link = linkiter
                nodedist(node,link) = sqrt((nodepos(node,1)-link_startpos(counter,1))^2+(nodepos(node,2)-link_startpos(counter,2))^2);
                nodedist(node,link+1) = sqrt((nodepos(node,1)-link_endpos(counter,1))^2+(nodepos(node,2)-link_endpos(counter,2))^2);
                counter = counter + 1;
            end
        end
        %Append each node to appropriate link and location
        for node = 1:size(nodeinds,1)
           tmp_dist = nodedist(node,:);
           [~,loc] = min(tmp_dist);
           loc_check = loc/2;
           if floor(loc_check)==loc_check %is loc_check an integer?
              link_add = loc_check; %Link to add node to
              link_start = 0; %Add it to the end
           else %Loc check isn't an integer
              link_add = ceil(loc_check); %Link to add node to
              link_start = 1; %Add it to the beginning
           end
           %Add to seglinkinds
           tmp_seglink = seglinkinds(link_add,:);
           tmp_seglink = tmp_seglink(~isnan(tmp_seglink));
           if link_start == 1
              tmp_seglink = [nonzeros(nodeinds(node,:))',tmp_seglink];
           else
              tmp_seglink = [tmp_seglink,nonzeros(nodeinds(node,:))'];
           end
           seglinkinds(link_add,1:length(tmp_seglink)) = tmp_seglink;
        end    
        seglinkinds(seglinkinds==0) = NaN; %Added to prevent anamalous 0's from popping up in seglinkinds (from matlab automatically extending matrix)
        %Done adding nodes

        if size(seglinktmp,1)>2
            avgpos = zeros(length(seglinktmp),2);
            for i=1:length(seglinktmp) %Calc subs for each point in each seg
                %     for i=1:1
                subtmp = zeros(length(seglinkinds(~isnan(seglinkinds(i,:)))),2);
                for j=1:size(subtmp,1)
                    [subtmp(j,1),subtmp(j,2)] = ind2sub(size(toplot.mask),seglinkinds(i,j));   %Get all subinds for the link
                end
                avgrow = mean(subtmp(:,1));
                avgcol = mean(subtmp(:,2));
                avgpos(i,1) = avgrow; %Save average position for each segment in the joined segment
                avgpos(i,2) = avgcol;
            end
            positions = size(avgpos,1);
            dist = zeros(positions,positions-1);
            for j=1:positions
                pos_vec = 1:positions;
                positer = pos_vec(~ismember(pos_vec,j));
                counter = 1;
                for i = positer
                    dist(j,counter) = sqrt((avgpos(i,2)-avgpos(j,2))^2 + (avgpos(i,1)-avgpos(j,1))^2);
                    counter = counter + 1;
                end
            end            
            [dist,dist_sortinds] = sort(dist,2); %Now sort the rows of the dist matrix in rows
            sum_dist = sum(dist,2); %Sum distances in each row (starting link)
            [~,distind] = max(sum_dist); %Get the starting segment (one with max sum of distances)
            seglinkinds_tmp = seglinkinds;
            seglinkinds_tmp(1,:) = seglinkinds(distind,:); %this is now our first link
            positer = pos_vec(~ismember(pos_vec,distind));
            for i = 1:size(dist_sortinds,2)
                seglinkinds_tmp(i+1,:) = seglinkinds(positer(dist_sortinds(distind,i)),:); %switch around rows based on sorting
            end
            seglinkinds = seglinkinds_tmp;
            %Now order segments in order of increasing distance from #1, then flip
%             if necessary. Then replace in vsl_graph_tmp.
        end
        
        %Now need to flip links if their indices aren't in line with prev links
        seglinkinds_nonan = struct();
        for i=1:size(seglinkinds,1)
            tmpinds1 = seglinkinds(i,:);
            tmpinds2 = tmpinds1(~isnan(tmpinds1));
            seglinkinds_nonan(i).inds = tmpinds2;
        end
        joinedlinktmp = seglinkinds_nonan(1).inds;
        for i=2:size(seglinkinds,1)
            joinedlinktmp = [joinedlinktmp,seglinkinds_nonan(i).inds];
        end
        
        linklength = zeros(size(seglinkinds,1),1);
        for i=1:size(seglinkinds,1)
            linklength(i) = sum(~isnan(seglinkinds(i,:)));
        end
        counter = 1;
        for i=1:length(joinedlinktmp)-1
            [row1,col1] = ind2sub(size(im_mask),joinedlinktmp(i));
            [row2,col2] = ind2sub(size(im_mask),joinedlinktmp(i+1));
            dxtmp(i) = sqrt((row2-row1)^2+(col2-col1)^2);
        end
        
        % FLIP LINKS STARTING NOW
        %For first link, check if first pix is closer to link 2 than last
        %pix (then need to flip it). Could also do this for other links
        %instead of while.
        nextlinkinds = seglinkinds(2,:);
        [link2row,link2col] = ind2sub(size(im_mask),nextlinkinds);
        [pix1r,pix1c] = ind2sub(size(im_mask),seglinkinds(1,1)); %First pix of first segment
        [pix2r,pix2c] = ind2sub(size(im_mask),seglinkinds(1,sum(~isnan(seglinkinds(1,:))))); %Last pix of first segment
        dist1 = sqrt((pix1r-mean(link2row(~isnan(link2row))))^2+(pix1c-mean(link2col(~isnan(link2col))))^2);
        dist2 = sqrt((pix2r-mean(link2row(~isnan(link2row))))^2+(pix2c-mean(link2col(~isnan(link2col))))^2);
        if dist1<dist2 %Then we need to flip FIRST segment
            seglinkinds_toflip = seglinkinds(1,:);
            seglinkinds_toflip = seglinkinds_toflip(~isnan(seglinkinds_toflip));
            seglinkinds(1,1:length(seglinkinds_toflip)) = flip(seglinkinds_toflip);
            seglinkinds_nonan = struct();
            for i=1:size(seglinkinds,1)
                tmpinds1 = seglinkinds(i,:);
                tmpinds2 = tmpinds1(~isnan(tmpinds1));
                seglinkinds_nonan(i).inds = tmpinds2;
            end
        end
        joinedlinktmp = seglinkinds_nonan(1).inds;
        for i=2:size(seglinkinds,1)
            joinedlinktmp = [joinedlinktmp,seglinkinds_nonan(i).inds];
        end
        
        
        distthresh = 10; %Pixels (5 was too small)
        if sum(dxtmp>distthresh)>0 %If we need to flip any links
            for prevlink = 1:size(seglinkinds,1)-1
                prevlinkinds = seglinkinds(prevlink,:); %Compare avg previous link location to first and last of the next link (this has NaNs)
                [link2row,link2col] = ind2sub(size(im_mask),prevlinkinds);
                [pix1r,pix1c] = ind2sub(size(im_mask),seglinkinds(prevlink+1,1)); %First pixel of next link
                [pix2r,pix2c] = ind2sub(size(im_mask),seglinkinds(prevlink+1,sum(~isnan(seglinkinds(prevlink+1,:))))); %Last pixel of the next link
                prevlinkinds_last = seglinkinds(prevlink,sum(~isnan(prevlinkinds)));
                [link2rowlast,link2collast] = ind2sub(size(im_mask),prevlinkinds_last);
                dist3 = sqrt((link2rowlast-pix1r)^2+(link2collast-pix1c)^2); %Last pixel of first seg to first pix of next seg
                dist4 = sqrt((link2rowlast-pix2r)^2+(link2collast-pix2c)^2); %Last pixel of first seg to last pix of next seg
                if dist3>dist4 %Then we need to flip the next segment (Might need to restrict this part more..) like && dist1 and dist2 are close to eachother.
                    seglinkinds_toflip = seglinkinds(prevlink+1,:);
                    seglinkinds_toflip = seglinkinds_toflip(~isnan(seglinkinds_toflip));
                    seglinkinds(prevlink+1,1:length(seglinkinds_toflip)) = flip(seglinkinds_toflip);
                    seglinkinds_nonan = struct();
                    for i=1:size(seglinkinds,1)
                        tmpinds1 = seglinkinds(i,:);
                        tmpinds2 = tmpinds1(~isnan(tmpinds1));
                        seglinkinds_nonan(i).inds = tmpinds2;
                    end                    
                end
            end
            joinedlinktmp = seglinkinds_nonan(1).inds;
            for i=2:size(seglinkinds,1)
                joinedlinktmp = [joinedlinktmp,seglinkinds_nonan(i).inds];
            end     
        end
        %Update vsl_graph_tmp entries: (cc_ind)
        vsl_graph_tmp.link.cc_ind{startlink,1} = joinedlinktmp;
        if mod(seglinkiter,10)==0
            seglinkiter
        end
    end

    %Finally, delete links that were incorporated into others above.
    todel = [];
    for i=1:size(vsl_graph_tmp.link.seglinks,1)
        nonzeroseglinks = nonzeros(vsl_graph_tmp.link.seglinks(i,:));
        todel = [todel;nonzeroseglinks(2:end)]; %Delete second-through-end seglinks from cc_ind
    end
    tmptodel = zeros(length(vsl_graph_tmp.link.num_voxel_per_cc),1);
    tmptodel(todel) = 1;
    
    link_lengths = zeros(size(vsl_graph_tmp.link.cc_ind,1),1);
    for i=1:length(vsl_graph_tmp.link.cc_ind)
        if size(vsl_graph_tmp.link.cc_ind{i,1},1)==1
            vsl_graph_tmp.link.cc_ind{i,1} = vsl_graph_tmp.link.cc_ind{i,1}';
        end
        link_lengths(i) = size(vsl_graph_tmp.link.cc_ind{i,1},1);
    end
    single_pix = link_lengths == 1;
    tmptodel = or(tmptodel,single_pix);

    links = vsl_graph_tmp.link.cc_ind; %Linear indices of all links
    links(tmptodel) = []; %Delete short links that only have one node
    vsl_graph_tmp.link.cc_ind = []; %Replace cc_inds after deleting 
    vsl_graph_tmp.link.cc_ind = links;
    %Update other vsl_graph_tmp.link entries
    vsl_graph_tmp.link.pos_ind = []; %Replace list of all link indices
    vsl_graph_tmp.link.pos_ind = cell2mat(vsl_graph_tmp.link.cc_ind(:)); %Concatenate remaining link indices
    cnltmp = vsl_graph_tmp.link.connected_node_label;
    cnltmp(tmptodel,:) = []; %Delete node labels for deleted links
    vsl_graph_tmp.link.connected_node_label = cnltmp;
    [nrows,~] = cellfun(@size,vsl_graph_tmp.link.cc_ind);
    vsl_graph_tmp.link.num_voxel_per_cc = []; %Delete link lengths for deleted links
    vsl_graph_tmp.link.num_voxel_per_cc = nrows; %Number of segments in each link after joining and deleting
    numnodetmp = vsl_graph_tmp.link.num_node;
    numnodetmp(logical(tmptodel)) = []; %Replace number of nodes for each remaining link
    vsl_graph_tmp.link.num_node = numnodetmp;

    vsl_graph_tmp.link.num_cc = sum(~tmptodel); %Remaining links in the graph
    vsl_graph_tmp.link.num_voxel = size(vsl_graph_tmp.link.pos_ind,1);
    vsl_graph_tmp.link.label = repelem(1:vsl_graph_tmp.link.num_cc, vsl_graph_tmp.link.num_voxel_per_cc)';
    vsl_graph_tmp.link.map_ind_2_label = sparse(vsl_graph_tmp.link.pos_ind, ...
        ones(vsl_graph_tmp.link.num_voxel,1), ...
        vsl_graph_tmp.link.label, ...
        vsl_graph_tmp.num.block_voxel,1);
    %Now all "link" fields are updated as necessary
    
    for i=1:length(vsl_graph_tmp.link.cc_ind)
       linklengths(i) = length(vsl_graph_tmp.link.cc_ind{i,1}); 
    end
    figure 
    ecdf(vsl_graph.link.num_voxel_per_cc/pix_mm)
    hold on
    ecdf(linklengths/pix_mm)
    xlabel('Length (mm)')
    ylabel('Probability')
    title('Cumulative Probability Distribution Seg Lengths, old (blue), modified (red)')
    cd(savefile)
    savefig('ECDFSegmentLengths.fig');
end   

vsl_graph = vsl_graph_tmp;
clearvars -except link_lengths wind analyzed_files stim prev_link_results_struct vsl_graph files im_mask toplot wave savefile filename animal file pix_mm segs_link neurowave neurotoplot is_near_midline_Q rim
%Delete segments with only one pixel
link_lengths = zeros(size(vsl_graph.link.cc_ind,1),1);
for i=1:length(vsl_graph.link.cc_ind)
    if size(vsl_graph.link.cc_ind{i,1},1)==1
        vsl_graph.link.cc_ind{i,1} = vsl_graph.link.cc_ind{i,1}';
    end
    link_lengths(i) = size(vsl_graph.link.cc_ind{i,1},1);
end
vsl_graph.link.cc_ind(link_lengths==1) = [];
for i=1:length(vsl_graph.link.cc_ind)
    link_lengths(i) = size(vsl_graph.link.cc_ind{i,1},1);
end  

%% Phase Analysis
%Get first segment for all vsl_graph.liink.cc_ind cells
link_vec = zeros(length(vsl_graph.link.cc_ind),1);
link_skelval = zeros(length(vsl_graph.link.cc_ind),1);
link_wave = zeros(size(wave,2),length(vsl_graph.link.cc_ind)); %Chronux wants data in time x trials
link_loc = zeros(length(vsl_graph.link.cc_ind),2);

nodemat = zeros(200,length(vsl_graph.link.cc_ind)); %We don't know how many segments there will be yet, just preallocate 200 as an upper limit

for cell = 1:length(link_vec)
    link_vec(cell) = vsl_graph.link.cc_ind{cell,1}(1); %Get first link pixel
end
testmask = zeros(size(im_mask));
testmask(toplot.mask_ind) = 1;
for cell = 1:length(link_vec)
    if ~testmask(link_vec(cell)) %If the link pixel isn't in im_mask, find the small segment it belongs to
        disp('testmask error'); cell
        
        surmapind = NaN(8,2);
        mapind = link_vec(cell);
        surmapind(1,1) = mapind + 1;
        surmapind(2,1) = mapind - 1;
        surmapind(3,1) = mapind + size(toplot.mask,1);
        surmapind(4,1) = mapind - size(toplot.mask,1);
        surmapind(5,1) = mapind + size(toplot.mask,1) + 1;
        surmapind(6,1) = mapind + size(toplot.mask,1) - 1;
        surmapind(7,1) = mapind - size(toplot.mask,1) + 1;
        surmapind(8,1) = mapind - size(toplot.mask,1) - 1;
        for j = 1:8
            if ~isempty(toplot.skel_label(toplot.mask_ind==surmapind(j,1)))
                surmapind(j,2) = toplot.skel_label(toplot.mask_ind==surmapind(j,1));
            end
        end
        link_skelval(cell) = mode(surmapind(:,2));
        link_wave(:,cell) = wave(link_skelval(cell),:)'; %Get wave for this segment
        %Get x,y locations of all pixels in segment
        [row,col] = ind2sub(size(im_mask),toplot.mask_ind(toplot.skel_label==link_skelval(cell)));
        avgrow = round(mean(row));
        avgcol = round(mean(col));
        link_loc(cell,1) = avgrow;
        link_loc(cell,2) = avgcol;
        
    else
        link_skelval(cell) = toplot.skel_label(toplot.mask_ind==link_vec(cell)); %This is which segment the link pixel is in.
        link_wave(:,cell) = wave(link_skelval(cell),:)'; %Get wave for each segment
        %Get x,y locations of all pixels in segment
        [row,col] = ind2sub(size(im_mask),toplot.mask_ind(toplot.skel_label==link_skelval(cell)));
        avgrow = round(mean(row));
        avgcol = round(mean(col));
        link_loc(cell,1) = avgrow;
        link_loc(cell,2) = avgcol;
    end
end
link_results_struct = struct();
link_results_struct(1).wave = link_wave;
link_results_struct(1).locs = link_loc;


%Get the rest of the segment's wave and locations:
next_segs = struct();
maxlinklength = max(link_lengths);
segment_linknum = ones(length(vsl_graph.link.cc_ind),1);
link_results_struct(1).linknum = segment_linknum;
link_results_struct(1).skel_lab = link_skelval;
%Need to find the next link pixel that isn't in the previous segment
allNaNs = 0;
segment = 2;
while allNaNs == 0
    
    link_vec2 = zeros(length(vsl_graph.link.cc_ind),1);
    link_skelval2 = zeros(length(vsl_graph.link.cc_ind),1);
    link_wave2 = zeros(size(wave,2),length(vsl_graph.link.cc_ind)); %Chronux wants data in time x trials
    link_loc2 = zeros(length(vsl_graph.link.cc_ind),2);
    oneseg_link = zeros(length(vsl_graph.link.cc_ind),1);
    link_notin_maskind = zeros(length(vsl_graph.link.cc_ind),1);
    segment_linknum = link_results_struct(segment-1).linknum;
    prev_segment_linknum = segment_linknum; 

    for cell = 1:length(vsl_graph.link.cc_ind) %Iterate over vessels
        if length(vsl_graph.link.cc_ind{cell,1})>=segment_linknum(cell)+1
            link_vec2(cell) = vsl_graph.link.cc_ind{cell,1}(segment_linknum(cell)+1); %Get second link pixel
            if isempty(toplot.skel_label(toplot.mask_ind==link_vec2(cell))) %If the link pixel isn't in mask_ind, find the small segment it belongs to
                link_notin_maskind(cell) = 1;
            else
                repeat=1;
                nextlink = 1;
                while repeat==1 %Go through cc_ind pixels 
                    if ~isempty(toplot.skel_label(toplot.mask_ind==link_vec2(cell)))
                        link_skelval2(cell) = toplot.skel_label(toplot.mask_ind==link_vec2(cell)); %This is which segment the link pixel is in.
                        if link_skelval2(cell) == link_results_struct(segment-1).skel_lab(cell) %If we're still in the previous seg
                            if link_lengths(cell) >= prev_segment_linknum(cell) + 1 + nextlink %Can't go for longer than the link length
                                link_vec2(cell) = vsl_graph.link.cc_ind{cell,1}(prev_segment_linknum(cell) + 1 + nextlink); %Go down the link
                                repeat = 1; %We still need to repeat the above section
                                nextlink = nextlink + 1; %Keep getting consecuative link pixels until you reach a new small skel_label segment
                            else
                                oneseg_link(cell) = 1;
                                repeat = 0;
                            end
                        else
                            repeat = 0; %If we reached a new skel_label segment, stop repeating
                        end
                    else
                        repeat = 0;
                    end
                end
                segment_linknum(cell) = segment_linknum(cell)+nextlink; %If this link skelval~=prev one then segment_linknum = segment_linknum + 1. Otherwise add more nec. links
                
                link_wave2(:,cell) = wave(link_skelval2(cell),:)'; %Get wave for each segment
                %Get x,y locations of all pixels in segment
                [row,col] = ind2sub(size(im_mask),toplot.mask_ind(toplot.skel_label==link_skelval2(cell))); %Where this skel label is on the map
                avgrow = round(mean(row));
                avgcol = round(mean(col));
                link_loc2(cell,1) = avgrow;
                link_loc2(cell,2) = avgcol;

            end
            if segment_linknum(cell) == prev_segment_linknum(cell)
                segment_linknum(cell) = prev_segment_linknum(cell) + 1;
            end
            
        else
            link_wave2(:,cell) = NaN;
            link_loc2(cell,:) = NaN;
            segment_linknum(cell) = NaN;
            link_results_struct(segment).wave = link_wave2;
            link_results_struct(segment).locs = link_loc2;
            link_results_struct(segment).linknum = segment_linknum;
        end
        
    end
    link_wave2(:,find(oneseg_link)) = NaN; %Can change how we handle links that only span one segment (also can't calc. phase grad for only two segments R^2 = 1)
    
    linkinds = find(link_notin_maskind); %For link pixels not in map_ind, sample surrounding pixels and assign a segment based on their most common segment
    for i=1:length(linkinds)
        surmapind = NaN(8,2);
        linkind = linkinds(i);
        mapind = link_vec2(linkind);
        surmapind(1,1) = mapind + 1;
        surmapind(2,1) = mapind - 1;
        surmapind(3,1) = mapind + size(toplot.mask,1);
        surmapind(4,1) = mapind - size(toplot.mask,1);
        surmapind(5,1) = mapind + size(toplot.mask,1) + 1;
        surmapind(6,1) = mapind + size(toplot.mask,1) - 1;
        surmapind(7,1) = mapind - size(toplot.mask,1) + 1;
        surmapind(8,1) = mapind - size(toplot.mask,1) - 1;
        for j = 1:8
            if ~isempty(toplot.skel_label(toplot.mask_ind==surmapind(j,1)))
                surmapind(j,2) = toplot.skel_label(toplot.mask_ind==surmapind(j,1));
            end
        end
        link_skelval2(linkind) = mode(surmapind(:,2));
        link_wave2(:,linkind) = wave(link_skelval2(linkind),:)'; %Get wave for each segment
        %Get x,y locations of all pixels in segment
        [row,col] = ind2sub(size(im_mask),toplot.mask_ind(toplot.skel_label==link_skelval2(linkind)));
        avgrow = round(mean(row));
        avgcol = round(mean(col));
        link_loc2(linkinds(i),1) = avgrow;
        link_loc2(linkinds(i),2) = avgcol;
    end

    link_results_struct(segment).wave = link_wave2;
    link_results_struct(segment).locs = link_loc2;
    link_results_struct(segment).linknum = segment_linknum;
    link_results_struct(segment).skel_lab = link_skelval2;
    
    test = link_results_struct(segment).wave(1,:);
    testnan = isnan(test);
    numnan = sum(testnan);
    if numnan == length(vsl_graph.link.cc_ind) %
        allNaNs = 1;
    end
    
    str = sprintf('Segment %d done',segment);
    disp(str)
    segment = segment + 1;
end

%Now we have wave files for all segments, organized by each link they belong to.
%% Begin phase calculations
tic;

link_results_struct(1).file = files(file).name;
link_results_struct(1).folder = files(file).folder;

test = zeros(length(link_results_struct),length(vsl_graph.link.cc_ind));
for i=1:length(link_results_struct)
    for j=1:length(vsl_graph.link.cc_ind)
        test(i,j) = link_results_struct(i).wave(1,j);  %First wave value from all segments, all links (check for NaNs)     
    end
end
segs_link = zeros(length(vsl_graph.link.cc_ind),1);
for j=1:length(vsl_graph.link.cc_ind)
    segments_vec = test(:,j);
    non_nanstmp = ~isnan(segments_vec);
    segs_link(j) = sum(non_nanstmp); %Number of segments in each link
end
%Initialize CHRONUX
addpath(genpath('C:\chronux_2_12'));
t = 1:size(wave,2);
params.Fs = toplot.rate; %use this for interpolated data
params.pad    =    2;
params.fpass = [0 1];
% params.err   = [2 .05];
params.err   = 0;
params.trialave = 1; %Average across trials
T = size(wave,2)/toplot.rate; %We want T*BW ~= 5
if T>400
    vesBW = 0.01;
else
    vesBW = 0.02;
end
params.tapers = [round(T*vesBW),round(2*T*vesBW-1)] %Max number of tapers without distortion
for j=1:length(vsl_graph.link.cc_ind)    
    wavevec = zeros(size(wave,2),segs_link(j));
    for i=1:segs_link(j)
       wavevec(:,i) = link_results_struct(i).wave(:,j);         
    end
    [S,f]=mtspectrumc(wavevec,params);
    link_results_struct(1).AvgSpec(:,j) = S; %Spectra using Chronux averaging across vessel segment spectra
    link_results_struct(1).f(:,j) = f;
    if mod(j,10)==0
        str = sprintf('VesselSpec %.00f Percent Done',j/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
end
toc
spec_vec = link_results_struct(1).AvgSpec;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
for i=1:size(spec_vec,2)
    plot(f,spec_vec(:,i))
    hold on
end
subplot(2,1,2)
plot(f,log10(mean(spec_vec,2))); grid on;
[wind,~] = ginput(2);
findf1 = round(wind(1),3);
findf2 = round(wind(2),3);
link_results_struct(1).wind = wind; %Save frequency window used to find vasomotor peaks
%We have the window we'll use, now just find peaks for each link.

rmpath(genpath('C:\chronux_2_12'));
for j=1:length(vsl_graph.link.cc_ind)
    if j==1
        f1 = find(round(f,3)==round(findf1,3),1,'first');
        if isempty(f1)
            f1 = find(round(f,3)==round((findf1-0.001),3),1,'first'); %Push along findf to reach an f element
        end
        f2 = find(round(f,3)==round(findf2,3),1,'last');
        if isempty(f2)
            f2 = find(round(f,3)==round((findf2+0.001),3),1,'last');
        end
    end
    [maxes,f_peaks] = findpeaks(link_results_struct(1).AvgSpec(f1:f2,j),link_results_struct(1).f(f1:f2,j));
    maxloc = find(maxes==max(maxes));
    f_peak = f_peaks(maxloc);
    link_results_struct(1).f_peak(j) = f_peak;
    if mod(j,10)==0
        str = sprintf('VesselSpec2 %.00f Percent Done',j/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
end
%% Analyze neural data
clearvars cell
link_results_struct_neu.inds = cell(100,length(vsl_graph.link.cc_ind));

neu_rad = 8; %Make this changable parameter - pixels to extend from vessel edge (imdilate)
buffer_rad = 2; %Pixels, start neural signal away from vessels
tol = 4; %Tolerance to be used for segmenting neu_mask pixels

for i=1:length(vsl_graph.link.cc_ind)  
    %Need to assign radius to each link pixel
    mask_distance_transform = bwdist(~im_mask,'euclidean'); %From Xiang's fun_graph_get_voxel_radius_from_mask
    vsl_graph.radius.cc_ind_r{i,1} = mask_distance_transform(vsl_graph.link.cc_ind{i,1}); %Store radius for every link pixel 
    vsl_graph.radius.cc_ind_medr(i) = median(vsl_graph.radius.cc_ind_r{i,1}); %Store link median radius
    %Create dilated vessel
    neumap = zeros(size(im_mask));
    skel_lab_tmp = zeros(length(link_results_struct),1);
    for ii = 1:length(link_results_struct)
       skel_lab_tmp(ii) = link_results_struct(ii).skel_lab(i); 
    end
    skel_lab_tmp = nonzeros(skel_lab_tmp); %Skeleton labels for the current link
    skel_toplot = ismember(toplot.skel_label,skel_lab_tmp);
    neumap(toplot.mask_ind(skel_toplot)) = 1;
    neumap = imfill(neumap); %Fill in any gaps within the vessel
    
    neumap_dist = imdilate(neumap,strel('disk',neu_rad));
    neu_mask = double(im_mask);
    buff_mask = imdilate(neu_mask,strel('disk',buffer_rad));
    mask_tot = buff_mask;
    %Need to delete parts of dilated vessel (vessel and edges)
    vsl_graph.radius.cc_ind{i,1} = find(neumap_dist.*~mask_tot); %Surrounding pixels for each link 
    if mod(i,10) == 0
        str = sprintf('Neuromap %.00f Percent Done',i/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
end

%Now we need to segment the neural cc_ind into the same #segments per link
%as the corresponding vessel
%% 
clearvars segs pixlength step segiter skel_tmp j i ang_vec run_tot link_results_struct_neu.inds inds_tot map
 link_results_struct_neu(1).inds = cell(200,length(vsl_graph.link.cc_ind));
for i=1:length(vsl_graph.link.cc_ind)
    segs = segs_link(i); %Segments in current cc_ind link    
    pixlength =  length(vsl_graph.link.cc_ind{i,1});%Pixels in the link
    if segs >= 3 && pixlength >= 5
    step = pixlength/segs;
    segiter = round(1:step:pixlength);
    %Make dividing line each pixlength/segs

    skel_tmp = vsl_graph.link.cc_ind{i,1};
    rad_tmp = vsl_graph.radius.cc_ind_r{i,1};
    for j=1:length(segiter)-1        
        if j==1
            linkpix_tmp = skel_tmp(1);
            radpix_tmp = double(rad_tmp(1));
            [endp1_row,endp1_col] = ind2sub(size(im_mask),linkpix_tmp);
            surrp1 = skel_tmp(2:5);
            [surrp1_rows,surrp1_cols] = ind2sub(size(im_mask),surrp1);
            ang = atan2((mean(surrp1_rows)-endp1_row),(mean(surrp1_cols)-endp1_col));
            delx = (radpix_tmp + neu_rad + tol)*sin(ang);
            dely = (radpix_tmp + neu_rad + tol)*cos(ang);
            xi = [endp1_col+delx,endp1_col-delx];
            yi = [endp1_row-dely,endp1_row+dely];
            delx2 = 2*(radpix_tmp + neu_rad + tol)*cos(ang);
            dely2 = 2*(radpix_tmp + neu_rad + tol)*sin(ang);
            xi(3:4) = [xi(2)-delx2,xi(1)-delx2];
            yi(3:4) = [yi(2)-dely2,yi(1)-dely2];
            neubox = poly2mask(xi,yi,size(im_mask,1),size(im_mask,2));
            inds_tot = (vsl_graph.radius.cc_ind{i,1}(ismember(vsl_graph.radius.cc_ind{i,1},find(neubox))))';
            run_tot = inds_tot;
            link_results_struct_neu.inds{j,i} = (vsl_graph.radius.cc_ind{i,1}(ismember(vsl_graph.radius.cc_ind{i,1},find(neubox))))';
        else
            skelpix = segiter(j);
            radpix_tmp = double(rad_tmp(j));
            linkpix_tmp = skel_tmp(skelpix);
            [endp1_row,endp1_col] = ind2sub(size(im_mask),linkpix_tmp);
            if j < length(segiter)-1
                if length(skel_tmp) >= skelpix + 3
                surrp1 = skel_tmp(skelpix:skelpix + 3);
                else
                   surrp1 = skel_tmp(skelpix:end); 
                end
            else
                surrp1 = skel_tmp(skelpix:end); 
            end
            [surrp1_rows,surrp1_cols] = ind2sub(size(im_mask),surrp1);
            ang = atan2((mean(surrp1_rows)-endp1_row),(mean(surrp1_cols)-endp1_col));
            delx = (radpix_tmp + neu_rad + tol)*sin(ang);
            dely = (radpix_tmp + neu_rad + tol)*cos(ang);
            xi = [endp1_col+delx,endp1_col-delx];
            yi = [endp1_row-dely,endp1_row+dely];
            delx2 = 2*(radpix_tmp + neu_rad + tol)*cos(ang);
            dely2 = 2*(radpix_tmp + neu_rad + tol)*sin(ang);
            xi(3:4) = [xi(2)-delx2,xi(1)-delx2];
            yi(3:4) = [yi(2)-dely2,yi(1)-dely2];
            if j < length(segiter)-1
                neubox = poly2mask(xi,yi,size(im_mask,1),size(im_mask,2));
                inds_tot = (vsl_graph.radius.cc_ind{i,1}(ismember(vsl_graph.radius.cc_ind{i,1},find(neubox))))';
                link_results_struct_neu.inds{j,i} = inds_tot(~ismember(inds_tot,run_tot));
                run_tot = [run_tot,inds_tot(~ismember(inds_tot,run_tot))];
            else %If this is the last line, get all remaining pixels
                inds_tot = vsl_graph.radius.cc_ind{i,1};
                link_results_struct_neu.inds{j,i} = inds_tot(~ismember(inds_tot,run_tot));
                run_tot = [run_tot,inds_tot(~ismember(inds_tot,run_tot))'];
            end
        end
    end
    if mod(i,10) == 0
        str = sprintf('Dividing neural segs %.00f Percent Done',i/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
    end
end
%Now refine the neural graph
%Get rid of disconnected elements and empty cells
for i = 1:length(vsl_graph.link.cc_ind)
    nonempty = ~cellfun(@isempty,link_results_struct_neu(1).inds(:,i));
    nsegs_tmp = link_results_struct_neu(1).inds(nonempty,i);
    for j=1:length(nsegs_tmp)
        map = zeros(size(im_mask));
        map(nsegs_tmp{j}) = 1;
        map = bwareaopen(map,2,4);
        nsegs_tmp{j} = find(map);
    end
    nonempty2 = ~cellfun(@isempty,nsegs_tmp);
    nsegs_tmp2 = nsegs_tmp(nonempty2);
    link_results_struct_neu(1).inds(1:sum(nonempty2),i) = nsegs_tmp2;
    endempty = cell(size(link_results_struct_neu(1).inds,1)-sum(nonempty2),1);
    link_results_struct_neu(1).inds(sum(nonempty2)+1:size(link_results_struct_neu.inds,1),i) = endempty;
    if mod(i,10) == 0
        str = sprintf('Cleaning neural segs %.00f Percent Done',i/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
end


map = zeros(size(im_mask)); %Plot neuron mask
for ii = 1:size(link_results_struct_neu(1).inds,2) %Update this for other masks
    for j=1:size(link_results_struct_neu(1).inds,1)
        map(link_results_struct_neu(1).inds{j,ii}) = j;
    end
end
figure
imagesc(map+double(im_mask)*-10)
caxis([-10 50]);
colormap jet
daspect([1,1,1]);
colorbar

%% Get wave time series for each neural segment.
n_segs_link = zeros(size(link_results_struct_neu(1).inds,2),1);
for i=1:length(n_segs_link) 
    counter = 0; 
    for j=1:size(link_results_struct_neu(1).inds)
       if ~isempty(link_results_struct_neu(1).inds{j,i})
          counter = counter + 1; 
       end
    end
    n_segs_link(i) = counter;
end

%Populate wave values
for link = 1:length(vsl_graph.link.cc_ind) %Iterate over vessels
    nonempty = ~cellfun(@isempty,link_results_struct_neu(1).inds(:,link));
    if link == 1
        for j=1:max(n_segs_link)
            link_results_struct_neu(j).wave = NaN(size(neurowave,2),size(link_results_struct_neu(1).inds,2));
        end
    end
    for j=1:max(n_segs_link) %Iterate over segments
        if nonempty(j)==1
            inds_tmp = cell2mat(link_results_struct_neu(1).inds(j,link));
            skel_inds = ismember(neurotoplot.mask_ind,inds_tmp);
            nwave_tmp = neurowave(neurotoplot.skel_label(skel_inds),:)';
            link_results_struct_neu(j).wave(:,link) = mean(nwave_tmp,2);
        else
            link_results_struct_neu(j).wave(:,link) = NaN;
        end
    end
    if mod(link,10) == 0
        str = sprintf('Getting neural waves %.00f Percent Done',link/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
end

%% Calculate phase differences and distance for the neural segments %%

tic
addpath(genpath('C:\chronux_2_12'));
t = 1:size(neurowave,2);
params.Fs = neurotoplot.rate; %use this for interpolated data
params.pad    =    2;
params.fpass = [0 1];
params.err   = [2 .05];
params.trialave = 1; %Try averaging across trials (check with manual averaging)
T = size(neurowave,2)/neurotoplot.rate;
if T>400 && T<800
    BW = 0.03;
elseif T<400
    BW = 0.06;
elseif T>800
    BW = 0.015;
end
params.tapers = [round(T*BW),round(2*T*BW-1)] %Max number of tapers without distortion
link_results_struct_neu(1).neuparams = params;
%Put together the arrays
[tapers1,pad1,FsStim1,fpass1,err1,trialave1]=getparams(params);
N1 = size(neurowave,2);
nfft1=max(2^(nextpow2(N1)+pad1),N1);
[f1,findx1]=getfgrid(FsStim1,nfft1,fpass1);
tapers1=dpsschk(tapers1,N1,FsStim1);

n_phi_mat = NaN(max(n_segs_link),size(link_results_struct_neu(1).wave,2));
n_dist_mat = NaN(max(n_segs_link),size(link_results_struct_neu(1).wave,2));
n_phi_mat_fv = NaN(max(n_segs_link),size(link_results_struct_neu(1).wave,2));
n_phi_mat_fs = NaN(max(n_segs_link),size(link_results_struct_neu(1).wave,2));
n_phi_mat_fs2 = NaN(max(n_segs_link),size(link_results_struct_neu(1).wave,2));
for i=1:max(n_segs_link)
    link_results_struct_neu(i).locs = NaN(size(link_results_struct_neu(1).wave,2),2);
    link_results_struct_neu(i).skel_locs = NaN(size(link_results_struct_neu(1).wave,2),2);
end

tpf_peak = toplot.f_peak;
vaso_findx_tpfv = max(find(round(f1,3)==round(tpf_peak(1),3))); %Global peak vasofrequency from beforepuff
vaso_findx_fs = max(find(round(f1,3)==round(tpf_peak(2),3))); %Stim frequency
vaso_findx_fs2 = max(find(round(f1,3)==round(tpf_peak(3),3))); %Stim 2nd harmonic
if isempty(vaso_findx_tpfv) || isempty(vaso_findx_fs) || isempty(vaso_findx_fs2)
    vaso_findx_tpfv = max(find(round(f1,2)==round(tpf_peak(1),2))); %Global peak vasofrequency from beforepuff
    vaso_findx_fs = max(find(round(f1,2)==round(tpf_peak(2),2))); %Stim frequency
    vaso_findx_fs2 = max(find(round(f1,2)==round(tpf_peak(3),2))); %Stim 2nd harmonic
end
link_results_struct(1).tpf_peak = tpf_peak; %this is [freq @ peak before-stim vasomotor power, stim freq, stim freq 2nd harmonic]

pix_wind = 10; %Tunable, this is to avoid jumping across parts of the vessels that are curved
for j=1:size(link_results_struct_neu(1).wave,2) %Iterate over links
    counter = 1;
    for seg = 1:n_segs_link(j)
        inds_tmp = cell2mat(link_results_struct_neu(1).inds(seg,j));
        [rtmp,ctmp] = ind2sub(size(im_mask),inds_tmp);
        avgr = mean(rtmp);
        avgc = mean(ctmp);
        link_results_struct_neu(seg).locs(j,1) = avgr;
        link_results_struct_neu(seg).locs(j,2) = avgc;
        
        if ~isnan(link_results_struct_neu(seg).locs(j,1)) && seg==1 %Make skel_locs; for first and last neural segs, use average position (extends past the skeleton).
            link_results_struct_neu(seg).skel_locs(j,:) = round(link_results_struct_neu(seg).locs(j,:));
        elseif ~isnan(link_results_struct_neu(seg).locs(j,1)) && seg>1 && seg<n_segs_link(j)
            if length(vsl_graph.link.cc_ind{j,1}) >= (counter + pix_wind) %If this isn't a small number of pixels
                skelinds = vsl_graph.link.cc_ind{j,1}(counter:counter + pix_wind);
            else
                skelinds = vsl_graph.link.cc_ind{j,1};
            end
            [rtmp,ctmp] = ind2sub(size(im_mask),skelinds);
            distvec = sqrt((rtmp-link_results_struct_neu(seg).locs(j,1)).^2+(ctmp-link_results_struct_neu(seg).locs(j,2)).^2);
            [~,I] = min(distvec);
            mindist_skelind = skelinds(I);
            %update counter
            counter = find(vsl_graph.link.cc_ind{j,1}==mindist_skelind); %Index in cc_ind of current skel link
            [link_results_struct_neu(seg).skel_locs(j,1),link_results_struct_neu(seg).skel_locs(j,2)] = ind2sub(size(im_mask),mindist_skelind);
        elseif ~isnan(link_results_struct_neu(seg).locs(j,1)) && seg == n_segs_link(j)
            link_results_struct_neu(seg).skel_locs(j,:) = round(link_results_struct_neu(seg).locs(j,:));
        end
    end
        
    vaso_findx = max(find(round(f1,3)==round(link_results_struct(1).f_peak(j),3))); %Vessel-wise vasofreq
    for i=1:n_segs_link(j) %Iterate over segments
        if i==1 %If we're on the first segment
            wave0 = link_results_struct_neu(i).wave(:,j);
            Jneu0 = mtfftc(wave0,tapers1,nfft1,FsStim1);
            Jneuplot0 = Jneu0(findx1,:,:);
            clear Jves0

            wave1 = link_results_struct_neu(i).wave(:,j);
            Jneu1 = mtfftc(wave1,tapers1,nfft1,FsStim1);
            Jneuplot1 = Jneu1(findx1,:,:);
            clear Jves1
            
            S12 = squeeze(mean(conj(Jneuplot1(vaso_findx,:)).*Jneuplot0(vaso_findx,:),2));
            n_phi_mat(i,j) = angle(S12);
            S12_fv = squeeze(mean(conj(Jneuplot1(vaso_findx_tpfv,:)).*Jneuplot0(vaso_findx_tpfv,:),2));
            n_phi_mat_fv(i,j) = angle(S12_fv);
            S12_fs = squeeze(mean(conj(Jneuplot1(vaso_findx_fs,:)).*Jneuplot0(vaso_findx_fs,:),2));
            n_phi_mat_fs(i,j) = angle(S12_fs);
            S12_fs2 = squeeze(mean(conj(Jneuplot1(vaso_findx_fs2,:)).*Jneuplot0(vaso_findx_fs2,:),2));
            n_phi_mat_fs2(i,j) = angle(S12_fs2);
            n_dist_mat(i,j) = 0;
        else
            wave1 = link_results_struct_neu(i).wave(:,j);
            Jneu1 = mtfftc(wave1,tapers1,nfft1,FsStim1);
            Jneuplot1 = Jneu1(findx1,:,:);
            clear Jves1
            
            S12=squeeze(mean(conj(Jneuplot1(vaso_findx,:)).*Jneuplot0(vaso_findx,:),2));
            n_phi_mat(i,j) = angle(S12);
            S12_fv = squeeze(mean(conj(Jneuplot1(vaso_findx_tpfv,:)).*Jneuplot0(vaso_findx_tpfv,:),2));
            n_phi_mat_fv(i,j) = angle(S12_fv);
            S12_fs = squeeze(mean(conj(Jneuplot1(vaso_findx_fs,:)).*Jneuplot0(vaso_findx_fs,:),2));
            n_phi_mat_fs(i,j) = angle(S12_fs);
            S12_fs2 = squeeze(mean(conj(Jneuplot1(vaso_findx_fs2,:)).*Jneuplot0(vaso_findx_fs2,:),2));
            n_phi_mat_fs2(i,j) = angle(S12_fs2);
            
            dx = (link_results_struct_neu(i).skel_locs(j,2)-link_results_struct_neu(i-1).skel_locs(j,2));
            dy = (link_results_struct_neu(i).skel_locs(j,1)-link_results_struct_neu(i-1).skel_locs(j,1));
            n_dist_mat(i,j) = sqrt(dx^2+dy^2)+n_dist_mat(i-1,j);
        end
    end  
    if mod(j,10)==0
        str = sprintf('Neural PhaseDiff %.00f Percent Done',j/length(vsl_graph.link.cc_ind) * 100);
        disp(str);
    end
end
toc


%% Calculate vessel phase differences %%
tic
addpath(genpath('C:\chronux_2_12'));
t = 1:size(wave,2);
params.Fs = toplot.rate; %use this for interpolated data
params.pad    =    2;
params.fpass = [0 1];
params.err   = [2 .05];
params.trialave = 1; %Try averaging across trials (check with manual averaging)
T = size(wave,2)/toplot.rate;
if T>400 && T<800
    BW = 0.03;
elseif T<400
    BW = 0.06;
elseif T>800
    BW = 0.015;
end
params.tapers = [T*BW,round(2*T*BW-1)]; %Max number of tapers without distortion
link_results_struct(1).params = params;
%Put together the arrays
[tapers1,pad1,FsStim1,fpass1,err1,trialave1]=getparams(params);
N1 = size(wave,2);
nfft1=max(2^(nextpow2(N1)+pad1),N1);
[f1,findx1]=getfgrid(FsStim1,nfft1,fpass1);
tapers1=dpsschk(tapers1,N1,FsStim1);

phi_mat = NaN(max(segs_link),length(vsl_graph.link.cc_ind));
dist_mat = NaN(max(segs_link),length(vsl_graph.link.cc_ind));
phi_mat_fv = NaN(max(segs_link),size(link_results_struct(1).wave,2));
phi_mat_fs = NaN(max(segs_link),size(link_results_struct(1).wave,2));
phi_mat_fs2 = NaN(max(segs_link),size(link_results_struct(1).wave,2));
node_mat = NaN(size(phi_mat));

for j=1:length(vsl_graph.link.cc_ind) %Iterate over links
    vaso_findx = max(find(round(f1,3)==round(link_results_struct(1).f_peak(j),3)));
    for i=1:segs_link(j) %Iterate over segments
        if i==1 %If we're on the first segment
            wave0 = link_results_struct(i).wave(:,j);
            Jves0 = mtfftc(wave0,tapers1,nfft1,FsStim1);
            Jvesplot0 = Jves0(findx1,:,:);
            clear Jves0
            
            wave1 = link_results_struct(i).wave(:,j);
            Jves1 = mtfftc(wave1,tapers1,nfft1,FsStim1);
            Jvesplot1 = Jves1(findx1,:,:);
            clear Jves1
            
            S12=squeeze(mean(conj(Jvesplot1(vaso_findx,:)).*Jvesplot0(vaso_findx,:),2));
            phi_mat(i,j) = angle(S12);
            S12_fv = squeeze(mean(conj(Jvesplot1(vaso_findx_tpfv,:)).*Jvesplot0(vaso_findx_tpfv,:),2));
            phi_mat_fv(i,j) = angle(S12_fv);
            S12_fs = squeeze(mean(conj(Jvesplot1(vaso_findx_fs,:)).*Jvesplot0(vaso_findx_fs,:),2));
            phi_mat_fs(i,j) = angle(S12_fs);
            S12_fs2 = squeeze(mean(conj(Jvesplot1(vaso_findx_fs2,:)).*Jvesplot0(vaso_findx_fs2,:),2));
            phi_mat_fs2(i,j) = angle(S12_fs2);
            dist_mat(i,j) = 0;

        else
            wave1 = link_results_struct(i).wave(:,j);
            Jves1 = mtfftc(wave1,tapers1,nfft1,FsStim1);
            Jvesplot1 = Jves1(findx1,:,:);
            clear Jves1
            
            S12=squeeze(mean(conj(Jvesplot1(vaso_findx,:)).*Jvesplot0(vaso_findx,:),2));
            phi_mat(i,j) = angle(S12);
            S12_fv = squeeze(mean(conj(Jvesplot1(vaso_findx_tpfv,:)).*Jvesplot0(vaso_findx_tpfv,:),2));
            phi_mat_fv(i,j) = angle(S12_fv);
            S12_fs = squeeze(mean(conj(Jvesplot1(vaso_findx_fs,:)).*Jvesplot0(vaso_findx_fs,:),2));
            phi_mat_fs(i,j) = angle(S12_fs);
            S12_fs2 = squeeze(mean(conj(Jvesplot1(vaso_findx_fs2,:)).*Jvesplot0(vaso_findx_fs2,:),2));
            phi_mat_fs2(i,j) = angle(S12_fs2);

            dx = (link_results_struct(i).locs(j,2)-link_results_struct(i-1).locs(j,2));
            dy = (link_results_struct(i).locs(j,1)-link_results_struct(i-1).locs(j,1));
            dist_mat(i,j) = sqrt(dx^2+dy^2)+dist_mat(i-1,j);
        end
    end  
    if mod(j,10)==0
        str = sprintf('PhaseGrad %.00f Percent Done',j/length(vsl_graph.link.cc_ind) * 100);
        disp(str)
    end
end
toc

dist_diff = diff(dist_mat);
max_diffs = max(dist_diff,[],'omitnan');
err_todel = max_diffs>(10*(pix_mm/90.4)); %Make sure link_results_struct locs and everything are updated. (apply to everything with 290 entries)

n_dist_diff = diff(n_dist_mat);
n_max_diffs = max(n_dist_diff,[],'omitnan');
n_err_todel = n_max_diffs>(20*pix_mm/115);

graph_errors = sum(err_todel) %This is how many segments weren't correlty joined/aligned
n_graph_errors = sum(n_err_todel)
link_results_struct(1).graph_errors = graph_errors;
link_results_struct_neu(1).n_graph_errors = n_graph_errors;

link_results_struct(1).phi_mat = phi_mat;
link_results_struct(1).phi_mat_fv = phi_mat_fv;
link_results_struct(1).phi_mat_fs = phi_mat_fs;
link_results_struct(1).phi_mat_fs2 = phi_mat_fs2;

link_results_struct(1).dist_mat = dist_mat;
link_results_struct(1).nodemat = nodemat;
link_results_struct_neu(1).n_phi_mat = n_phi_mat;
link_results_struct_neu(1).n_phi_mat_fv = n_phi_mat_fv;
link_results_struct_neu(1).n_phi_mat_fs = n_phi_mat_fs;
link_results_struct_neu(1).n_phi_mat_fs2 = n_phi_mat_fs2;
link_results_struct_neu(1).n_dist_mat = n_dist_mat;
link_results_struct_neu(1).f_peak = link_results_struct(1).f_peak; %Copy vessel vasomotion frequencies NEED TO DELETE ERRORS FROM BOTH NOW

err_todel = or(err_todel, n_err_todel);
total_graph_errors = sum(err_todel)
%Delete elements of link_results_struct whose vessel graph has an error
for i=1:length(link_results_struct)
   link_results_struct(i).wave(:,err_todel) = [];
   link_results_struct(i).locs(err_todel,:) = [];
   link_results_struct(i).linknum(err_todel) = [];
   link_results_struct(i).skel_lab(err_todel) = [];
end
for i=1:length(link_results_struct_neu)
    link_results_struct_neu(i).wave(:,err_todel) = [];
    link_results_struct_neu(i).locs(err_todel,:) = [];
%     link_results_struct_neu(i).wave_filt(:,err_todel) = [];
    link_results_struct_neu(i).skel_locs(err_todel,:) = [];
end
link_results_struct(1).AvgSpec(:,err_todel) = [];
link_results_struct(1).f(:,err_todel) = [];
link_results_struct(1).f_peak(err_todel) = [];
link_results_struct_neu(1).f_peak(err_todel) = [];
link_results_struct(1).phi_mat(:,err_todel) = [];
link_results_struct(1).phi_mat_fv(:,err_todel) = [];
link_results_struct(1).phi_mat_fs(:,err_todel) = [];
link_results_struct(1).phi_mat_fs2(:,err_todel) = [];
link_results_struct(1).dist_mat(:,err_todel) = [];
vsl_graph.link.cc_ind(err_todel,:) = [];
segs_link(err_todel) =  [];
n_segs_link(err_todel) = [];
link_results_struct_neu(1).n_segs_link = n_segs_link;
link_results_struct(1).segs_link = segs_link;

link_results_struct_neu(1).inds(:,err_todel) = [];
vsl_graph.radius.cc_ind_r(err_todel) = [];
vsl_graph.radius.cc_ind_medr(err_todel) = [];
vsl_graph.radius.cc_ind(err_todel) = [];
link_results_struct_neu(1).n_phi_mat(:,err_todel) = [];
link_results_struct_neu(1).n_phi_mat_fv(:,err_todel) = [];
link_results_struct_neu(1).n_phi_mat_fs(:,err_todel) = [];
link_results_struct_neu(1).n_phi_mat_fs2(:,err_todel) = [];
link_results_struct_neu(1).n_dist_mat(:,err_todel) = [];

nodemat(:,err_todel) = [];
link_results_struct(1).nodemat(:,err_todel) = [];
%Delete corresponding vessels from nodemat

%% Calculate Vessel Phase gradient (~1 second)
tic
RSC = 1;
[phasevec] = MakePhasevec(vsl_graph,link_results_struct(1).phi_mat,link_results_struct(1).dist_mat,RSC,pix_mm,link_results_struct);
[phasevec_fv] = MakePhasevec(vsl_graph,link_results_struct(1).phi_mat_fv,link_results_struct(1).dist_mat,RSC,pix_mm,link_results_struct);
[phasevec_fs] = MakePhasevec(vsl_graph,link_results_struct(1).phi_mat_fs,link_results_struct(1).dist_mat,RSC,pix_mm,link_results_struct);
[phasevec_fs2] = MakePhasevec(vsl_graph,link_results_struct(1).phi_mat_fs2,link_results_struct(1).dist_mat,RSC,pix_mm,link_results_struct);

link_results_struct(1).phasevec = phasevec;
link_results_struct(1).phasevec_fv = phasevec_fv;
link_results_struct(1).phasevec_fs = phasevec_fs;
link_results_struct(1).phasevec_fs2 = phasevec_fs2;
link_results_struct(1).vsl_graph = vsl_graph;

%Now neurons
[n_phasevec] = MakePhasevec(vsl_graph,link_results_struct_neu(1).n_phi_mat,link_results_struct_neu(1).n_dist_mat,RSC,pix_mm,link_results_struct);
[n_phasevec_fv] = MakePhasevec(vsl_graph,link_results_struct_neu(1).n_phi_mat_fv,link_results_struct_neu(1).n_dist_mat,RSC,pix_mm,link_results_struct);
[n_phasevec_fs] = MakePhasevec(vsl_graph,link_results_struct_neu(1).n_phi_mat_fs,link_results_struct_neu(1).n_dist_mat,RSC,pix_mm,link_results_struct);
[n_phasevec_fs2] = MakePhasevec(vsl_graph,link_results_struct_neu(1).n_phi_mat_fs2,link_results_struct_neu(1).n_dist_mat,RSC,pix_mm,link_results_struct);

link_results_struct_neu(1).n_phasevec = n_phasevec;
link_results_struct_neu(1).n_phasevec_fv = n_phasevec_fv;
link_results_struct_neu(1).n_phasevec_fs = n_phasevec_fs;
link_results_struct_neu(1).n_phasevec_fs2 = n_phasevec_fs2;
link_results_struct_neu(1).n_vsl_graph = vsl_graph;
toc
%% Get the node locations from locs field
nodemat = NaN(200,length(vsl_graph.link.cc_ind));
for i=1:length(vsl_graph.link.cc_ind)
    nodevec = zeros(segs_link(i),1);
    for j=1:segs_link(i)
        segind = sub2ind(size(im_mask),link_results_struct(j).locs(i,1),link_results_struct(j).locs(i,2));
        %Find which segment this belongs to
        find_mask_ind = find(toplot.mask_ind==segind);
        if ~isempty(find_mask_ind)
            skel(j) = toplot.skel_label(find_mask_ind);

            skelinds = find(toplot.skel_label == skel(j));
            mapinds = toplot.mask_ind(skelinds);
            membercheck = ismember(mapinds,vsl_graph.node.pos_ind);
            if any(membercheck)
                nodevec(j) = 1; %This (ith) vessel's (segment'th) has a node
            end
        end
    end
    nodemat(1:length(nodevec),i) = nodevec;
end

link_results_struct(1).nodemat = nodemat;
%% Perform t-test to exclude phase gradients with correlation coefficient not statistically different from 0. Reference Biostatistical Analysis p383
% This t-test is recalcluated in other analyses at difference significance
% levels. This test is equivalent to testing if the regression slope is
% statistically different from 0.

alpha = 0.05; %Set significance level
p = 1-alpha;
p2 = 1-alpha/2;
t_test_mat = zeros(length(vsl_graph.link.cc_ind),2);
for i=1:length(vsl_graph.link.cc_ind)
    r = link_results_struct(1).phasevec(i,2);
    n = segs_link(i);
    df = n-2;
    if df > 0
        SE = sqrt((1-r^2)/(n-2));
        t = r/SE;
        tcrit = icdf('T',p2,df);
        t_test_mat(i,1) = t;
        t_test_mat(i,2) = tcrit;
    else
        t_test_mat(i,1) = NaN;
        t_test_mat(i,2) = NaN;
    end
end
link_results_struct(1).t_test_mat = t_test_mat;

neu_t_test_mat = zeros(length(vsl_graph.link.cc_ind),2);
for i=1:length(vsl_graph.link.cc_ind)
    r = link_results_struct_neu(1).n_phasevec(i,2);
    n = n_segs_link(i);
    df = n-2;
    if df > 0
        SE = sqrt((1-r^2)/(n-2));
        t = r/SE;
        tcrit = icdf('T',p2,df);
        neu_t_test_mat(i,1) = t;
        neu_t_test_mat(i,2) = tcrit;
    else
        neu_t_test_mat(i,1) = NaN;
        neu_t_test_mat(i,2) = NaN;
    end
end
link_results_struct_neu(1).neu_t_test_mat = neu_t_test_mat;

%% Calculate vessel-wise v-n coherence and save

segs_linkv = link_results_struct(1).segs_link;
segs_linkn = link_results_struct_neu(1).n_segs_link;
TrialT = size(link_results_struct(1).wave,1);
TrialTn = size(link_results_struct_neu(1).wave,1);
Tdiff = TrialT-TrialTn;
%Calc mean wave signals for each link
for link = 1:length(segs_linkv)
    if Tdiff == 0

        wavemat = zeros(TrialT,segs_linkv(link));
        nwavemat = zeros(size(link_results_struct_neu(1).wave,1),segs_linkn(link));
        for seg = 1:segs_linkv(link)
            wavemat(:,seg) = link_results_struct(seg).wave(:,link);
        end
        for segn = 1:segs_linkn(link)
            nwavemat(:,segn) = link_results_struct_neu(segn).wave(:,link);
        end
        nwavenan = isnan(nwavemat(1,:));
        nwavemat(:,nwavenan) = [];
        avgwave = mean(wavemat,2);
        link_results_struct(1).avgwave(:,link) = avgwave;
        link_results_struct_neu(1).avgwave(:,link) = mean(nwavemat,2);
    elseif Tdiff > 0 %Vessels longer than neurons
        wavemat = zeros(TrialT-Tdiff,segs_linkv(link));
        nwavemat = zeros(TrialTn,segs_linkn(link));
        for seg = 1:segs_linkv(link)
            wavemat(:,seg) = link_results_struct(seg).wave(1:end-Tdiff,link);
        end
        for segn = 1:segs_linkn(link)
            nwavemat(:,segn) = link_results_struct_neu(segn).wave(:,link);
        end
        nwavenan = isnan(nwavemat(1,:));
        nwavemat(:,nwavenan) = [];
        avgwave = mean(wavemat,2);
        link_results_struct(1).avgwave(:,link) = avgwave;
        link_results_struct_neu(1).avgwave(:,link) = mean(nwavemat,2);
        N1 = max(size(link_results_struct_neu(1).avgwave(:,link)));
    elseif Tdiff < 0
        wavemat = zeros(TrialT,segs_linkv(link));
        nwavemat = zeros(TrialTn+Tdiff,segs_linkn(link));
        for seg = 1:segs_linkv(link)
            wavemat(:,seg) = link_results_struct(seg).wave(:,link);
        end
        for segn = 1:segs_linkn(link)
            nwavemat(:,segn) = link_results_struct_neu(segn).wave(1:end+Tdiff,link);
        end
        nwavenan = isnan(nwavemat(1,:));
        nwavemat(:,nwavenan) = [];
        avgwave = mean(wavemat,2);
        link_results_struct(1).avgwave(:,link) = avgwave;
        link_results_struct_neu(1).avgwave(:,link) = mean(nwavemat,2);
        N1 = max(size(link_results_struct_neu(1).avgwave(:,link)));
    end
end

%Calculate coherence
addpath(genpath('C:\chronux_2_12'));
t = 1:(TrialT-abs(Tdiff));
params.Fs = neurotoplot.rate;
params.pad   = 2;
params.fpass = [0 1];
params.err   = [2 .05];
params.trialave = 0;
if Tdiff<0
    T = TrialT/neurotoplot.rate;
elseif Tdiff>0
    T = (TrialT-Tdiff)/neurotoplot.rate;
else
    T = TrialT/neurotoplot.rate; 
end
if T>400 && T<800
    BW = 0.03;
elseif T<400
    BW = 0.06;
elseif T>800
    BW = 0.015;
end
params.tapers = [round(T*BW),round(2*T*BW-1)] %Max number of tapers without distortion
link_results_struct(1).cohparams = params;
%Put together the arrays
[NW,pad1,FsStim1,fpass1,err1,trialave1]=getparams(params);
nfft1=max(2^(nextpow2(N1)+pad1),N1);
[f1,findx1]=getfgrid(FsStim1,nfft1,fpass1);
tapers1=dpsschk(NW,N1,FsStim1);
link_results_struct(1).freqforcoh = f1;

C12_f = NaN(length(f1),length(segs_linkv));
% [C,phi,S12,S1,S2,f,ConfC,phistd,Cerr]=coherencyc(wavemat,nwavemat,params); %Same result as below code

for link = 1:length(segs_linkv) %Iterate over links in the trial
    vaso_freq = link_results_struct(1).f_peak(link); %Vessel vasomotor frequency
    if T>200
        vaso_findx = max(find(round(f1,3)==round(vaso_freq,3)));  %Find vasomotion frequency
    else
        vaso_findx = max(find(round(f1,2)==round(vaso_freq,2)));  %Find vasomotion frequency
    end
    
    Jves0 = mtfftc(link_results_struct(1).avgwave(:,link),tapers1,nfft1,FsStim1);
    Jvesplot0 = Jves0(findx1,:,:);
    if ~isnan(link_results_struct_neu(1).avgwave(1,link))
        Jneu0 = mtfftc(link_results_struct_neu(1).avgwave(:,link),tapers1,nfft1,FsStim1);
        Jneuplot0 = Jneu0(findx1,:,:);
        
        S12=squeeze(mean(conj(Jvesplot0(vaso_findx,:)).*Jneuplot0(vaso_findx,:),2));
        S1=squeeze(mean(conj(Jvesplot0(vaso_findx,:)).*Jvesplot0(vaso_findx,:),2));
        S2=squeeze(mean(conj(Jneuplot0(vaso_findx,:)).*Jneuplot0(vaso_findx,:),2));
        C12=S12./sqrt(S1.*S2); %Coherence at vasomotion frequency
        phi=angle(C12);
        
        S12=squeeze(mean(conj(Jvesplot0(:,:)).*Jneuplot0(:,:),2));
        S1=squeeze(mean(conj(Jvesplot0(:,:)).*Jvesplot0(:,:),2));
        S2=squeeze(mean(conj(Jneuplot0(:,:)).*Jneuplot0(:,:),2));
        C12_f(:,link)=S12./sqrt(S1.*S2); %Coherence at all freqs
        phi_f=angle(C12_f(:,link));
        
        link_results_struct(1).vncoh(link,1) = C12; %Save single results at vaso freq
        link_results_struct(1).vnphi(link,1) = phi;
%         link_results_struct(1).vncoh_f(1:length(C12_f),link) = C12_f;
%         link_results_struct(1).vnphi_f(1:length(phi_f),link) = phi_f;
%         %Dont save these, make files large/long to load
        
        %Calculate correlation too
%         [C,lags] = xcorr(link_results_struct(1).avgwave(:,link),link_results_struct_neu(1).avgwave(:,link),'normalized');
%         link_results_struct(1).Corr(1:length(C),link) = C;
%         link_results_struct(1).lags(1:length(lags),link) = lags;
    else
        link_results_struct(1).vncoh(link,1) = NaN; %Save single results at vaso freq
        link_results_struct(1).vnphi(link,1) = NaN;
%         link_results_struct(1).vncoh_f(:,link) = NaN;
%         link_results_struct(1).vnphi_f(:,link) = NaN;
%         link_results_struct(1).Corr(1:(length(link_results_struct(1).avgwave(:,link))-1),link) = NaN;
%         link_results_struct(1).lags(1:(length(link_results_struct(1).avgwave(:,link))-1),link) = NaN;
    end
end 
avgcoh_vasofreq = mean(abs(link_results_struct(1).vncoh),'omitnan')%Single number
avgcoh_f = mean(abs(C12_f),2,'omitnan'); %Coherence as a fxn of f
link_results_struct(1).C12_f = C12_f;
link_results_struct(1).avgcoh_f = avgcoh_f;
% link_results_struct(1).avgvncoh = avgcoh; %Mean coherence at the
% vasomotion frequency Can just calculate if needed, don't save.

%% Calculate directionality using locations of segments
%Function to calculate direction of each vessel in link_results_struct
if ~strcmp(animal,'JD221024F2') && ~strcmp(animal,'TB200217F2')%This is an open window animal & signel hemisphere animal
[ves_ang] = VesselDirection(link_results_struct,is_near_midline_Q,rim);
link_results_struct(1).ves_ang = ves_ang;
end

%% Calculate phase and k as a function of time using Hilbert Transform
disp('Bandpass and Hilbert K')

PassF = 0.18;
StopF = 0.2;
if size(wave,2)/toplot.rate > 200
lowp = designfilt('lowpassfir','FilterOrder',200,'PassbandFrequency',PassF,'StopbandFrequency',StopF,'SampleRate',toplot.rate);
filtwave = filtfilt(lowp,wave'); %Wave needs to be in time x space here
else
   lowp = designfilt('lowpassfir','FilterOrder',100,'PassbandFrequency',PassF,'StopbandFrequency',StopF,'SampleRate',toplot.rate);
filtwave = filtfilt(lowp,wave'); %Wave needs to be in time x space here 
end
wavehphase = angle(hilbert(filtwave)); %Phase(t) x space

% For each vessel, calculate k(t) from phase(t)
tic;
dist_mat = link_results_struct(1).dist_mat;
hkmat = zeros(size(wavehphase,1),length(segs_linkv));
hRmat = zeros(size(wavehphase,1),length(segs_linkv));
htmat = zeros(size(wavehphase,1),length(segs_linkv));
htcrit = zeros(1,length(segs_linkv));
for i=1:length(segs_linkv) %Takes 17 seconds on CPU
    if segs_linkv(i) > 1
        numsegs = segs_linkv(i);
        skel_tmp = zeros(numsegs,1);
        for j=1:segs_linkv(i) %Iterate over segments in current vessel to get skel labels
            skel_tmp(j) = link_results_struct(j).skel_lab(i);
        end
        hphasemat = wavehphase(:,skel_tmp);
        distvec = dist_mat(:,i);
        distvec = distvec(~isnan(distvec));
        ktmp = zeros(size(wavehphase,1),1);
        Rtmp = zeros(size(wavehphase,1),1);
        for j = 1:size(wavehphase,1) %Iterate over time
            phivec = hphasemat(j,:);
            %Unwrap phase to correct for -pi -> pi discont
            tphase = mean(phivec);
            ophase = phivec - tphase;
            ophase = mod(ophase + pi, 2*pi) - pi;

            fit2 = polyfit(distvec,ophase,1);
            ktmp(j) = fit2(1);
            R = corrcoef(distvec,ophase);
            Rtmp(j) = R(1,2);
        end
        df = numsegs-2;
        if df > 0
            SE = ((1-Rtmp.^2)./(n-2)).^(1/2);
            ttmp = Rtmp./SE;
            tcrit = icdf('T',p2,df);
        else
            ttmp = NaN;
            tcrit = NaN;
        end
        %Save results
        hkmat(:,i) = ktmp;
        hRmat(:,i) = Rtmp;
        htmat(:,i) = ttmp;
        htcrit(i) = tcrit;
    end
end

HilbertK = struct();
HilbertK(1).kmat = hkmat;
HilbertK(1).Rmat = hRmat;
HilbertK(1).ktmat = htmat;
HilbertK(1).tcrit  = htcrit;
HilbertK(1).Lowp = [PassF,StopF];

link_results_struct(1).HilbertK = HilbertK;
toc



%% Analyze/filter/save results
%MAKE FUNCTION fun_saveandplot.m
for i=1:length(vsl_graph.link.cc_ind)
   max_vec(i) = max(link_results_struct(1).dist_mat(:,i)); 
   n_max_vec(i) = max(link_results_struct_neu(1).n_dist_mat(:,i));
end
link_results_struct(1).link_lengths_mm = max_vec/pix_mm;
link_results_struct_neu(1).n_link_lengths_mm = n_max_vec/pix_mm;

%Save link_results files before deleting anything
if isfield(toplot,f_peak)
    link_results_struct(1).meanvasopeak = toplot.f_peak(1);
    link_results_struct(1).fwind = [findf1,findf2];
    link_results_struct_neu(1).meanvasopeak = toplot.f_peak(1);
    link_results_struct_neu(1).fwind = [findf1,findf2];
end

save([filename,'results.mat'],'link_results_struct');
save([filename,'neural_results.mat'],'link_results_struct_neu')

plotQ = 0; %1 to plot trial summary figures
cd(savefile)
fun_plotandsave(plotQ,link_results_struct,link_results_struct_neu,t_test_mat,neu_t_test_mat,vsl_graph,max_vec,n_max_vec,pix_mm,phasevec,n_phasevec,f,animal,wind)

file
end

