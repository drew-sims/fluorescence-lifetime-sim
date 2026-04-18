function [t, x, h, y] = fluorescence_emmision(tau, TL, TS, Tmax)
    % Inputs:
    % tau  : Fluorescence lifetime (ns)
    % TL   : Laser repetition period (ns)
    % TS   : Sampling period (ns)
    % Tmax : Total simulation time (ns)

    % 1. Create time vector
    t = 0:TS:Tmax;
    
    x = zeros(size(t));
    pulse_indices = 1:round(TL/TS):length(t);
    x(pulse_indices) = 1;

    h = exp(-t / tau); 
    
    y_full = conv(x, h);
    y = y_full(1:length(t));
end