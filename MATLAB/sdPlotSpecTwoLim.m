function sdPlotSpecTwoLim(f, g, gsd, llim, ulim, fignum)
%Quick function to plot a spectrum with upper and lower SD spectra.
%
%sdPlotSpec(f, g, gsd, llim, ulim, fignum)
%
%INPUTS
%   f   - Frequency axis
%   g   - Spectral signal
%   gsd - StD of spectrum
%   llim, ulim - Plots in frequency window of ulim - llim
%   fignum -   Desired figure number window.

ig = find((f < ulim) - (f < llim));
figure(fignum);semilogy(f(ig), abs(g(ig))); hold; semilogy(f(ig), abs(g(ig)) + abs(gsd(ig)), 'red');
semilogy(f(ig), abs(g(ig)) - abs(gsd(ig)), 'green');
xlabel('F (Hz'); ylabel('abs power'); title('Spectra with SD spectra');
legend('Signal Spectra', 'Signal + 1 SD', 'Signal - 1 SD');