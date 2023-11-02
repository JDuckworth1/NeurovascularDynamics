% fun_get_foldername.m

function [folder_name] = fun_get_foldername(filename)

if contains(filename,'2019')
    folder_name =string(extractAfter(filename,'2019_'));
elseif contains(filename,'2020')
    folder_name =string(extractAfter(filename,'2020_'));
elseif contains(filename,'2021')
    folder_name =string(extractAfter(filename,'2021_'));
elseif contains(filename,'2022')
    folder_name =string(extractAfter(filename,'2022_'));
else
    folder_name = [];
    disp('Trial year not found')
end

end



