function [fs, Ps] = FFT_PMTs(PMTchan, dt);
%Gives FFT of a PMT channel.
%Inputs:
%   PMTchan - All data from a single PMT over a period of time
%   dt - time element ( 1/ f)
%Outputs:
%   fs - Frequencies
%   Ps - (Complex) transform values
xc = xcorr(PMTchan - mean(PMTchan));
sxc = size(xc)
w = gausswin(sxc(2), 40);
%[fs,Ps] = spec(xc .* w.', dt);
[fs,Ps] = spec(xc, dt);