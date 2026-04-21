function hf = get_fluorescence_irf(t, tau)
    % hf(t) = exp(-t/tau) * u(t)
    hf = exp(-t/tau);
    hf(t < 0) = 0;
end