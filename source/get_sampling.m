function [y_n, t_sampled] = get_sampling(t, y_t, Ts)
    
    t_sampled = 0:Ts:max(t);
    
    y_n = interp1(t, y_t, t_sampled, 'nearest', 0);

end