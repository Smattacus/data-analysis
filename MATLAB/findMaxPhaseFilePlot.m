function [phase1, phase2] = findMaxPhaseFilePlot(fn, Ta, Fs, Fc, plots)
%Function which calculates and plots the maximum phase for the square wave.
%This assumes that you are using the synchronization box, and generates a
%square wave programmatically rather than using inverse transforms.
%
%It will display plots of the phase for the user to verify that the phase
%is correct and the same for both lasers.
%
%[phase1, phase2] = findMaxPhase(filename, N, Fs, Fc, plots)
%
d = h5read(fn, '/PMT_DATA_8BIT');
s1 = sum(d(7:9,:));
s2 = sum(d(23:25,:));
phases = linspace(0, 2 * pi * Fc * Ta, Ta * Fs);
p = linspace(0, 2 * pi , 100);
PMT1 = zeros(1, 100);
PMT2 = zeros(1, 100);
for i=1:100
    PMT1(i) = sum((square(phases + p(i)) + 1) .* s1) / 2;
    PMT2(i) = sum((square(phases + p(i)) + 1) .* s2) / 2;
end
% if plots == true
%     figure(1); plot(p, PMT1); xlabel('Phase (rad)'); ylabel('a.u.'); title('PMT 1 Phase plot');
%     figure(2); plot(p, PMT2); xlabel('Phase (rad)'); ylabel('a.u.'); title('PMT 2 Phase plot');
% end
im1 = find(PMT1 == max(PMT1));
phase1 = p(im1);
im2 = find(PMT2 == max(PMT2));
phase2 = p(im2);
if plots==true
    h1 = figure(1); set(h1, 'visible', 'off');
    plot(p, PMT1); title(sprintf('PMT #1, \nSquare Wave Corr for file %s', fn));
    xlabel(' Displacement (radians)');
    ylabel('Correlation');
    xcfn = sprintf('PhasePNG/%s_PHASE_PMT1.png', strrep(fn, '.h5', ''));
    legend(sprintf('Max at %f rad', phase1));
    print(h1, '-dpng', xcfn);
    h2 = figure(2); set(h2, 'visible', 'off');
    plot(p, PMT2); title(sprintf('PMT #2,\nSquare Wave Corr for file %s', fn));
    xlabel('Displacement (2 pi / 12)');
    ylabel('Correlation');
    legend(sprintf('Max = %f', phase2));
    xcfn = sprintf('PhasePNG/%s_PHASE_PMT2.png', strrep(fn, '.h5', ''));
    print(h2, '-dpng', xcfn);
end