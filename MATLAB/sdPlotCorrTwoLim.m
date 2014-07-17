function sdPlotCorrTwoLim(t, corr, corrsd, llim, ulim, fignum)
%Quick function to plot a spectrum with upper and lower SD spectra.
%
%It is HIGHLY RECOMMENDED to close the desired figure beforehand.
%
%sdPlotSpecTwoLim(f, g, gsd, llim, ulim, fignum)
%
%INPUTS
%   t   - time
%   corr - Correlation signal
%   gsd - StD of spectrum
%   llim, ulim - Shows window given by ulim - llim
%   fignum -   Desired figure number window.
figure(fignum);
it = find((t < ulim) - (t < llim));
plot(t(it), corr(it)); 
hold;
plot(t(it), corr(it) + corrsd(it), 'red');
plot(t(it), corr(it) - corrsd(it), 'green');
xlabel('Time (s)'); ylabel('XCorr value');
title('Correlation function with standard deviations');
legend('Correlation Function', 'Corr + StD', 'Corr - StD');