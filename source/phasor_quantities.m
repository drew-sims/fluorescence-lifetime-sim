function [s, g, tau] = phasor_quantities(TL, t, h_f)

    wL = 2*pi/TL;
    
    s_f = h_f*sin(wL);
    g_f = h_f*cos(wL);
    
    format long
    denom = trapz(t, h_f);
    s = trapz(t, s_f) / denom;
    g = trapz(t, g_f) / denom;
    
    tau = s/(wL*g);
end