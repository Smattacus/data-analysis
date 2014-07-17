function plotSVDPCSpec(filename, tc, freq, max_graphs)
%Quick function to plot the spectra of the principal components of the SVD
%of a typical langmuir probe array data file.
%
%plotSVDPCSPec(filename, tc)
%
% INPUTS
%   filename - Input filename folllowing Feng's lanmuir probe array data
%   format.
%   tc - Correlation length for windowing. Set tc = 0 not to window.
%   freq - Sampling frequency of the data. This is USUALLY 55.5 kHz.
%   max_graphs - Highest mode # to graph. Set to 8 to display all of them.
if nargin < 4
    display('Number of Desired Modes not Entered! Outputting all 8...');
    max_graphs = 8;
end
ac = [];
[T, XC] = SvdTimeCorr(filename);
N = size(T(:, 1), 1);
if tc ~= 0
    t = (-(N-1):(N-1)) / freq;
    win = exp(-(t/tc).^2/2).';
else
    display('Warning: Spectra are not windowed.');
    win = 1;
end
close all
for i=1:max_graphs
   temp = xcorr(T(:,i) - mean(T(:,i)), 'unbiased');
   [f, g] = spec(temp .* win, 1/55e3);
   ig = find((f < 20000) - (f < 0));
   figure(i); semilogy(f(ig), abs(g(ig)));
   title(sprintf('Normalized Power Spectrum for Mode #%d', i));
   xlabel('Frequency');
   ylabel('Power');
end