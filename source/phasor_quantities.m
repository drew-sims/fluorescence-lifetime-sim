function [s, g, tau_approx] = phasor_quantities(TL, t, h_f)
    omega = 2 * pi / TL;   % fundamental frequency (rad/ns)

    denom = sum(h_f);

    % If pixel is empty, return zeros
    if denom == 0
        s = 0; g = 0; tau_approx = 0;
        return;
    end

    g = sum(h_f .* cos(omega * t)) / denom;
    s = sum(h_f .* sin(omega * t)) / denom;

    % Recover lifetime from phasor coordinates
    tau_approx = s / (omega * g);   % phase lifetime (ns)
end