function x = get_excitation(t, TL)
    % x(t) = sum of delta(t - n*TL)
    % We approximate the delta by setting the sample closest to n*TL to 1/Ts
    x = zeros(size(t));
    pulse_times = 0:TL:max(t);
    for pt = pulse_times
        [~, idx] = min(abs(t - pt));
        x(idx) = 1;
    end
end