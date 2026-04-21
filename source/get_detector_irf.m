function hd = get_detector_irf(t, t0, sigma_t)
    % hd(t) = exp(-(t-t0)^2 / sigma_t^2) * u(t)
    hd = exp(-(t - t0).^2 / sigma_t^2);
    hd(t < 0) = 0; % Causal constraint
    % Normalize area to 1 to preserve signal energy
    hd = hd / sum(hd); 
end