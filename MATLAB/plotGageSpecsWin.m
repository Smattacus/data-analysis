function [f, ga, gb] = plotGageSpecsWin(filename, fig1, fig2, win)
%Quick function which plots the windowed autocorrelation of the two gage
%channels from the GAGE scope from an hdf5 file.
%
%[f, ga, gb] = plotGageSpecsWin(filename, fig1, fig2, win)
%
%win = The frequency below which to display. Plots from 0 - win Hz.
a = h5read(filename, '/GAGE_CHAN_A');
b = h5read(filename, '/GAGE_CHAN_B');
xa = xcorr(double(a - mean(a)), 'unbiased');
xb = xcorr(double(b - mean(b)), 'unbiased');
N = size(a, 1);
t = (-(N-1):(N-1)).'/1e5;
[f, ga] = spec(xa .* exp(-(t/0.01).^2/2), 1e-5);
[f, gb] = spec(xb .* exp(-(t/0.01).^2/2), 1e-5);
ig = find((f < win) - (f < 0));
figure(fig1); semilogy(f(ig), abs(ga(ig))); 
title(strrep(sprintf('Spectrum of GAGE CHANNEL A\nFile = %s', filename), '_', '-'));
xlabel('Frequency'); ylabel('Abs Power');
figure(fig2); semilogy(f(ig), abs(gb(ig))); 
title(strrep(sprintf('Spectrum of GAGE CHANNEL B\nFile = %s', filename), '_', '-'));
xlabel('Frequency'); ylabel('Abs Power');