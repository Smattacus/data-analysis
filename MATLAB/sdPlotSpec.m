function sdPlotSpec(f, g, gsd, lim, fignum)
%Quick function to plot a spectrum with upper and lower SD spectra.
%
%sdPlotSpec(f, g, gsd, lim, fignum)
%
%INPUTS
%   f   - Frequency axis
%   g   - Spectral signal
%   gsd - StD of spectrum
%   lim - absolute limit of frequency to plot up to.
%   fignum -   Desired figure number window.

ig = find(abs(f) < lim);
figure(fignum);semilogy(f(ig), abs(g(ig))); hold; semilogy(f(ig), abs(g(ig)) + abs(gsd(ig)), 'red');
semilogy(f(ig), abs(g(ig)) - abs(gsd(ig)), 'green');
xlabel('F (Hz'); ylabel('abs power'); title('Spectra with SD spectra');
legend('Signal Spectra', 'Signal + 1 SD', 'Signal - 1 SD');