function [z_ave, z_int, p_errs_ave, p_errs_int] = plotPhasesWin(f, specs, errs, ...
    ulim, llim, plotulim, plotllim, specindex, plotabscissa, xstring)
%Takes the data given and makes 3 plots:
% 1) Plot showing the area under the curve that is integrated / averaged.
% 2) The averaged phases with errorbars.
% 3) The integrated phases with errorbars.
%
% [z_ave, z_int, p_errs_ave, p_errs_int] = plotPhasesWin(f, specs, errs, ...
%     ulim, llim, plotulim, plotllim, specindex, plotabscissa, xstring)
%
[z_ave, p_errs_ave] = AvePhase(f, specs, errs, ulim, llim);
[z_int, p_errs_int] = IntPhase(f, specs, errs, ulim, llim);
close all;
figure(1);
%Plot area used.
ifwin = find((f < ulim) - (f < llim));
fi = find((f < plotulim) - (f < plotllim));
plot(f(fi), abs(specs(specindex, fi)));
hold; area(f(ifwin), abs(specs(specindex, ifwin)));
xlabel('F (Hz'); ylabel('|P|');
title('Area used to calculate phases');
%Plot averaged phases with errorbars.
figure(2);
errorbar(plotabscissa, angle(z_ave), p_errs_ave);
xlabel(xstring);
ylabel('Phases (rad)'); title('Averaged Phase plot');
%Plot integrated phases with errorbars.
figure(3);
errorbar(plotabscissa, angle(z_int), p_errs_int);
xlabel(xstring);
ylabel('Phases');
title('Integrated Phases');