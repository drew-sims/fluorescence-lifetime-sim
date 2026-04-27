function plot_phasor(TL, t, H)
% H: [nx, ny, nt] intensity cube (512x512x98)
% Each pixel's decay is H(i,j,:)

    [nx, ny, nt] = size(H);

    H = double(H);

    % Reshape to [npix x nt] so phasor_quantities can be vectorized
    npix = nx * ny;
    H_flat = reshape(H, [npix, nt]);   % [npix x nt]

    % Preallocate
    g_all       = zeros(npix, 1);
    s_all       = zeros(npix, 1);
    tau_all     = zeros(npix, 1);
    total_counts = sum(H_flat, 2);     % photon count per pixel
    
    figure();
    % Iterate over all pixels
    for k = 1:npix
        h_f = H_flat(k, :);           % 1D decay for this pixel
        [s_all(k), g_all(k), tau_all(k)] = phasor_quantities(TL, t, h_f);
        if mod(k,100) == 0
            plot(t, h_f);
            title("Fluorescence Lifetime");
            ylabel("Bin Counts");
            xlabel("Time (ns)");
        end
    end

    % Mask low-count pixels (avoids noisy phasors from background)
    min_counts = 10;
    mask = total_counts >= min_counts;

    g_plot   = g_all(mask);
    s_plot   = s_all(mask);
    tau_plot = tau_all(mask);

    % Plot
    figure('Color', 'white');
    hold on;

    % Universal semicircle
    theta = linspace(0, pi, 500);
    plot(0.5 + 0.5*cos(theta), 0.5*sin(theta), ...
         'k-', 'LineWidth', 1.5, 'DisplayName', 'Universal semicircle');

    % Scatter: each point is one pixel, colored by tau
    scatter(g_plot, s_plot, 6, tau_plot, 'filled', ...
            'MarkerFaceAlpha', 0.4);

    colormap(parula);
    cb = colorbar;
    cb.Label.String = '\tau (ns, recovered)';
    cb.Label.FontSize = 12;
    clim([prctile(tau_plot, 5), prctile(tau_plot, 95)]);  % robust color range

    xlabel('g', 'FontSize', 13);
    ylabel('s', 'FontSize', 13);
    title('Phasor Plot — all pixels', 'FontSize', 13);
    legend('Location', 'northeast');
    axis equal; xlim([0 1]); ylim([0 0.6]);
    grid on; box on;
    set(gca, 'FontSize', 11);
end