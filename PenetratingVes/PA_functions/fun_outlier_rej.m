
function [stats] = fun_outlier_rej(stats,rate)

    delD_Thresh = 11; %um/s
    delD_T_Frame = delD_Thresh/rate;
    if ~isfield(stats,'RadondEq_Outl')
        d_tmp = stats.RadondEq;
        d_diff = diff(d_tmp);
%         figure
%         plot(abs(d_diff),'k')
%         yline(delD_T_Frame);

        find_diff = find(abs(d_diff) > delD_T_Frame) + 1;
        %Check for consecuative indices. This handles 1 or 2 consecuative
        %jumps. Remaining cases need to be handled manually. 
        find_diff_check = diff(find_diff);
        for i = 1:length(find_diff_check)
            if find_diff_check(i) == 1
                vq = interp1([1,4],[d_tmp(find_diff(i)-1),d_tmp(find_diff(i+1)+1)],[1,2,3,4]);
                d_tmp([find_diff(i),find_diff(i+1)]) = vq([2,3]);
            else
                d_tmp(find_diff(i)) = mean([d_tmp(find_diff(i)-1),d_tmp(find_diff(i)+1)]);
            end
        end
%         figure
%         plot(d_tmp); title('shallow')

        stats.RadondEq_Outl = d_tmp;
        allstats = stats;
        allstats.rate = rate;

        save([animal,'_PA',PA,'_',depth1,'_allstats.mat'],'allstats');
    end
