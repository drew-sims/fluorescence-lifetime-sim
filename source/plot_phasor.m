function plot_phasor(t, h_f)
    
    % We put this in a for loop for each pixel, with separate h_f functions
    [s, g, tau] = phasor_quantities(12.5,t,h_f);
    
    % Plot (only for a single point, should be repurposed)
    figure('Color', 'white');
    scatter(g, s, 120, tau, 'filled');   % [tau] is explicitly a 1x1 double vector
    colormap(parula);
    cb = colorbar;
    clim([tau*0.9, tau*1.1]);
    cb.Label.String = '\tau (recovered)';
    cb.Label.FontSize = 12;
    xlabel('g',  'FontSize', 13);
    ylabel('s',  'FontSize', 13);
    title('Phasor Point (g, s) colored by \tau', 'FontSize', 13);
    grid on; box on;
    set(gca, 'FontSize', 11);
end