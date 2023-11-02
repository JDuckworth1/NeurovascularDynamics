%PAThresholdRadon.m
%Perform Radon transform and thresholding in Radon space using the process
%developed by Drew et al.

%Uses gpuArray for median filtering

function [dEq,CircVec,imgSeq,vesselPerim,time] = PAThresholdRadon(thresholdFWHM,thresholdInv,imgSeqin,rect,rate,current_cell,pix_um)

%Set up constant parameters
connectivity = 8;
CircVec = NaN(size(imgSeqin,3),1);
radonRVec = NaN(size(imgSeqin,3),1);
Outlier_Vec = NaN(size(imgSeqin,3),1);
RadonMedFilt = 5;

%% Crop image sequence
disp('Cropping and median filtering images')
tic;
for i=1:size(imgSeqin,3)
    if i==1
        tmp_im = imgSeqin(:,:,i);
        tmp_im = gpuArray(tmp_im);
        tmp_im_filt = medfilt2(tmp_im,[RadonMedFilt,RadonMedFilt]); %GPU medfilt2 must be odd between 3 and 15

        tmp_im_crop = imcrop(tmp_im_filt,rect);
        tmp_im_mean = tmp_im_crop - mean(tmp_im_crop(:));
        imgSeq_tmp = NaN(size(tmp_im_mean,1),size(tmp_im_mean,2),size(imgSeqin,3));
        tmp_im_mean = gather(tmp_im_mean);
        imgSeq_tmp(:,:,i) = tmp_im_mean; %Mean subtracted frame
    else
        tmp_im = imgSeqin(:,:,i);
        tmp_im = gpuArray(tmp_im);
        tmp_im_filt = medfilt2(tmp_im,[RadonMedFilt,RadonMedFilt]);

        tmp_im_crop = imcrop(tmp_im_filt,rect);
        tmp_im_mean = tmp_im_crop - mean(tmp_im_crop(:));
        tmp_im_mean = gather(tmp_im_mean);
        imgSeq_tmp(:,:,i) = tmp_im_mean;
    end
end
toc
imgSeq = imgSeq_tmp; clearvars imgSeq_tmp;
disp('Done cropping and median filtering images')

%% Set up parameters
[nRows, nCols, nFrames] = size(imgSeq);
areaPixels = zeros(nFrames, 1);
thetaRange = 0:179;
nAngles = length(thetaRange);
nRadii = 2*ceil(norm(size(imgSeq(:, :, 1)) - ...
    floor((size(imgSeq(:, :, 1)) - 1)/2) - 1)) + 3;
imgNorm = zeros(nRadii, nAngles, nFrames);
idxEdges = zeros(2, nAngles, nFrames);
imgInv = zeros(size(imgSeq));
vesselMask = false(size(imgInv));
vesselPerim = zeros(size(imgInv));

frameRate = rate;
frameTime = 1/frameRate;
if current_cell == 1
    t0 = (0.25*frameTime); %Average time for frame 1
else
    t0 = 3*(0.25*frameTime); %Average time for frame 2
end
time = (t0:frameTime:(nFrames*frameTime)-(frameTime-t0))';

filterType = 'Hamming';
outputSize = max([nRows, nCols]);
argsIRadon = {thetaRange, filterType, outputSize};
doResizeInv = false;
colsToUseInv = 1:nCols;
rowsToUseInv = 1:nRows;
if nRows ~= nCols
    doResizeInv = true;
    ctrOrig = floor(([nRows, nCols]+1)/2);
    ctrInv = floor((outputSize+1)/2);
    adjPixels = ctrInv - ctrOrig;
    hasRoundingProblem = all(adjPixels == 0);
    if hasRoundingProblem
        adjPixels = abs([nRows, nCols] - outputSize);
    end
    dimToAdj = find(adjPixels);
    if dimToAdj == 1
        rowsToUseInv = (1:nRows) + adjPixels(dimToAdj);
    else
        colsToUseInv = (1:nCols) + adjPixels(dimToAdj);
    end
end

%% Loop through all the frames
disp('Performing Radon transforms')
tic;
parfor iFrame = 1:nFrames
    % Calculate the radon transform of the frame
    imgTrans = radon(imgSeq(:, :, iFrame), thetaRange);
    [nRowsT, ~] = size(imgTrans);

    % Normalize each column (angle) of the radon transform
    imgNorm(:, :, iFrame) = (imgTrans - ...
        repmat(min(imgTrans, [], 1), [nRowsT, 1]))./ ...
        repmat(max(imgTrans, [], 1) - min(imgTrans, [], 1), ...
        [nRowsT, 1]);

    % Threshold the image, fill in any holes, and extract out
    % only the largest contiguous area
    imgThresh = imgNorm(:, :, iFrame) >= thresholdFWHM;
    imgThresh = imfill(imgThresh, 8, 'holes');
    ccT = bwconncomp(imgThresh, 4);
    ssT = regionprops(ccT);
    [~, maxIdx] = max([ssT(:).Area]);
    imgThresh = (labelmatrix(ccT) == maxIdx);

    for jAngle = 1:nAngles
        % Find the 'FWHM' edges of the transformed image
        idxEdgesTmp = [...
            find(imgThresh(:, jAngle), 1, 'first'), ...
            find(imgThresh(:, jAngle), 1, 'last')];

        % Manually threshold the transformed image using the
        % edges defined by the 'FWHM'
        imgThreshRowTmp = false(nRadii, 1);
        idxToUse = idxEdgesTmp(1) : idxEdgesTmp(2);
        imgThreshRowTmp(idxToUse) = true;
        imgThresh(:, jAngle) = imgThreshRowTmp;
        idxEdges(:, jAngle, iFrame) = idxEdgesTmp;
    end

    % Invert the thresholded radon-transformed image, adjust
    % the size if necessary, and then normalize it
    imgInvTemp = iradon(imgThresh, argsIRadon{:}); 
    if doResizeInv
        imgInvTemp = imgInvTemp(rowsToUseInv, colsToUseInv);
    end
    imgInv(:,:,iFrame) = imgInvTemp./max(imgInvTemp(:));

    % Threshold the inverted image, and fill in any holes in
    % the thresholded, inverted image
    imgInvThresh = imgInv(:,:,iFrame) > thresholdInv;
    imgInvThresh = imfill(imgInvThresh, 'holes');
    imgInvThresh = bwmorph(imgInvThresh,'majority'); %Added to remove hanging edges

    % Calculate the area of the largest contiguous region
    cc = bwconncomp(imgInvThresh, connectivity);
    ss = regionprops(cc,'Circularity','Area','BoundingBox');
    [areaPixels(iFrame, 1), maxIdx] = max([ss(:).Area]);

    % Create a binary image showing only the largest area
    % identified above
    lm = labelmatrix(cc);
    vesselMask(:, :, iFrame) = (lm == maxIdx);
    vesselPerim(:, :, iFrame) = bwperim((lm == maxIdx)); %Perimeter pixels for plotting
    CircVec(iFrame) = ss(maxIdx).Circularity; %Use this to check for good fit & threshold values.

end
toc
areaUm2 = areaPixels.*((1/pix_um)^2); %Square um
dEq = (4*areaUm2/pi).^(1/2);
disp('Done performing Radon transforms')


end





