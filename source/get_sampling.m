function [y_n, t_sampled] = get_sampling(t, y_t, Ts)
    % t: High-resolution time vector (simulation grid)
    % y_t: The continuous-time signal (d(t) * g(t))
    % Ts: Target sampling period (e.g., 0.2 ns)

    % Calculate the simulation step size
    dt = t(2) - t(1);
    
    % Calculate how many simulation steps are in one Ts
    % We use round() to handle any floating point floating point noise
    samples_per_Ts = round(Ts / dt);
    
    % Pick every N-th sample starting from the first one
    % This is the discrete-time equivalent of multiplying by p(t)
    indices = 1:samples_per_Ts:length(y_t);
    
    y_n = y_t(indices);
    t_sampled = t(indices);
end