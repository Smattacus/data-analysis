function [complex_means, phase_errors] = AvePhase(f, specs, errs, ulim, llim)
% This function returns the average phase over a certain window
% given by win. It then returns the error according to the erorr
% propagation equation, assuming that the errors in both Re(z)
% and Im(z) are independent and normally distributed.
%
% [complex_means, phase_errors] = AvePhase(f, specs, errs, ulim, llim)
%
%Inputs:
% f         - Frequency array for spectra.
% specs     - [Nspecs x SizeSpec] Array of spectra to find phases for.
% errs      - Array of reduced errors for the spectra.
% ulim      - Lower limit to cut off the spectrum.
% llim      - Upper limit to cut off the spectrum.
%
%Use ulim and llim to pick the range in which the averages are calculated.
%So, for example, 
%
%AvePhase(specs, errs, 1000, 1500)
%
%Would return the average phase and error bars for each spectrum for
%the interval 1 kHz to 1.5 kHz
Ns = size(specs, 1);
ifwin = find((f < ulim) - (f < llim));
complex_means = zeros(Ns, 1);
phase_errrors = zeros(Ns, 1);
for i=1:Ns
    complex_means(i) = mean(specs(i,ifwin));
    phase_errors(i) = sqrt(sum(0 * specs(i, ifwin) + (errs(i) ./ ...
    abs(specs(i, ifwin)) / size(ifwin,2)).^2));
end
