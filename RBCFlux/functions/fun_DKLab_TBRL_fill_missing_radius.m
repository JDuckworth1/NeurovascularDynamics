function r = fun_DKLab_TBRL_fill_missing_radius(r)


is_valid_Q = isfinite(r);
is_valid_t_Q = all(is_valid_Q, 1);
invalid_ind = find(~is_valid_t_Q);
if ~isempty(invalid_ind)
    fprintf('Number of invalid time points: %d\n', numel(invalid_ind));
    [num_trace, num_points] = size(r);
    % Replace data in all the columns by interpolation 
    x = 1 : num_points;
    x_v = x(is_valid_t_Q);
    x_q = x(invalid_ind);
    for iter_trace = 1 : num_trace
        ts_values = r(iter_trace, :);
        tmp_itp = griddedInterpolant(x_v, ts_values(is_valid_t_Q), 'linear', 'nearest');
        r(iter_trace, invalid_ind) = tmp_itp(x_q);    
    end
else
    fprintf('All the time points have valid value\n');
end
end
