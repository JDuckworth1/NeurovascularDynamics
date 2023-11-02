% fun_calc_ttestmat.m


function [t_test_mat] = fun_calc_ttestmat(pvcomb,alpha)

p = 1-alpha;
p2 = 1-alpha/2;
t_test_mat = zeros(length(pvcomb),2);
for i=1:length(pvcomb)
    r = pvcomb(i,2); %R value
    n = pvcomb(i,11); %Segs/link
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

end