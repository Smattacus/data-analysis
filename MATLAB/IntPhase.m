function [complex_ints, phase_errors] = IntPhase(f, specs, errs, ulim, llim)
% This function returns the integrated phase over a certain window
% given by win. It then returns the error according to the erorr
% propagation equation, assuming that the errors in both Re(z)
% and Im(z) are independent and normally distributed.
%
% [phases, phase_errors] = AvePhase(f, specs, errs, ulim, llim)
%
%Inputs:
% f         - Frequency array for spectra.
% specs     - [Nspecs x SizeSpec] Array of spectra to find phases for.
% errs      - Array of reduced errors for the spectra.
% ulim      - Lower limit to cut off the spectrum.
% llim      - Upper limit to cut off the spectrum.
%
%Use ulim and llim to pick the range in which the sums are calculated.
%So, for example, 
%
%IntPhase(specs, errs, 1000, 1500)
%
%Would return the integrated phase and error bars for each spectrum for
%the interval 1 kHz to 1.5 kHz
Ns = size(specs, 1);
ifwin = find((f < ulim) - (f < llim));
complex_ints = zeros(Ns, 1);
phase_errrors = zeros(Ns, 1);
for i=1:Ns
    complex_ints(i) = sum(specs(i,ifwin));
    phase_errors(i) = sqrt(sum(0 * specs(i, ifwin) + (errs(i) ./ ...
    abs(specs(i, ifwin))).^2));
end