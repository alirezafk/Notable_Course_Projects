function [h, h_range] = cross_corr(g, gi, c, ci, dt)
    h = dt*xcorr(g, c);
    h_length = length(h) - abs(length(g)-length(c));
    h = h(1:h_length);
    h_range = (0:(h_length-1))*dt + gi + ci;
    