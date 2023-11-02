%MakePhasevec.m

function [phasevec] = MakePhasevec(vsl_graph,phi_mat,dist_mat,RSC,pix_mm,link_results_struct)

phasevec = NaN(length(vsl_graph.link.cc_ind),8);
for i=1:length(vsl_graph.link.cc_ind)
    if sum(~isnan(phi_mat(:,i)))>2 %Need at least 3 points to fit the line
        phivec = zeros(sum(~isnan(phi_mat(:,i))),1);
        distvec = zeros(sum(~isnan(phi_mat(:,i))),1);
        for j=1:sum(~isnan(phi_mat(:,i)))
            phivec(j) = phi_mat(j,i);
            distvec(j) = dist_mat(j,i);
        end
        if RSC==1
            % RANSAC
            sampleSize = 2; % minimum size of sample from data required by polyfit to fit a model.
            maxDistance = 0.15; % max allowable distance for inliers. Right now this is arbitrary, max distance per point away from fit line.
            fitLineFcn = @(points) polyfit(distvec,phivec,1); % fit function using polyfit
            evalLineFcn = @(model, points) sum((phivec - polyval(model, distvec)).^2,2); % distance evaluation function, was sum(,2) which gave vector of values.
            [~, inlierIdx,status] = ransac([distvec,phivec],fitLineFcn,evalLineFcn, ...
                sampleSize,maxDistance);
            modelInliers = polyfit(distvec(inlierIdx),phivec(inlierIdx),1);
%             inlierPtsx = distvec(inlierIdx); %We don't need this 
%             inlierPtsy = phivec(inlierIdx);
%             x = [min(inlierPtsx) max(inlierPtsx)];
%             y = modelInliers(1)*x + modelInliers(2);
        end

        fit2 = polyfit(distvec,phivec,1);
        yhat = polyval(fit2,distvec);
        R = corrcoef(distvec,phivec); %Alternative is to use R2 formula, is the same as corrcoef for linear fits
        SSE = sum((phivec-yhat).^2);
%         ybar = sum(phivec)/numel(phivec); %Could use this if calculating R2 directly
        xbar = sum(distvec)/numel(distvec);
        Sxx = sum((distvec-xbar).^2);
        sigma = sqrt(numel(distvec)*SSE/((numel(distvec)-2)*Sxx));

        rad_pix = fit2(1);
        phasevec(i,1) = rad_pix*pix_mm; %rad_pix*pix_mm = k rad/mm
        phasevec(i,2) = R(1,2);
        phasevec(i,3) = link_results_struct(1).f_peak(i);
        phasevec(i,4) = sigma; %STD of slope (phase gradient)
        phasevec(i,5) = distvec(end); %Length of vessel segment
        phasevec(i,6) = i; %index of vsl_graph.link.cc_ind
%         if RSC == 1
        if status == 0 %No errors in the ransac process (see documentation), otherwise keep as NaN
%             if length(phivec) > 3
                phasevec(i,7) = modelInliers(1)*pix_mm; %k after RANSAC outlier rejection
%             end
        end
        phasevec(i,8) = sqrt(SSE/(numel(distvec)-2)); %Alternative measure of goodness of fit to R2
    end

end

end