function sdPlotCorr(t, corr, corrsd, lim, fignum)
%Quick function to plot a spectrum with upper and lower SD spectra.
%
%It is HIGHLY RECOMMENDED to close the desired figure beforehand.
%
%sdPlotSpec(f, g, gsd, lim, fignum)
%
%INPUTS
%   t   - time
%   corr - Correlation signal
%   gsd - StD of spectrum
%   lim - absolute limit of time to plot up to.
%   fignum -   Desired figure number window.
figure(fignum);
it = find(abs(t) < lim);
plot(t(it), corr(it)); 
hold;
plot(t(it), corr(it) + corrsd(it), 'red');
plot(t(it), corr(it) - corrsd(it), 'green');
xlabel('Time (s)'); ylabel('XCorr value');
title('Correlation function with standard deviations');
legend('Correlation Function', 'Corr + StD', 'Corr - StD');