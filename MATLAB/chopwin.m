function [chop_on chop_off] = chopwin(A, low, high, f);
%Creates two logical matrices that give chop_on and chop_off for a 
%
%  [chop_on chop_off] = chopwin(A, low, high, f)
%
%given dataset.
%Inputs:
%   A - Dataset
%   low - low limit on passband, in Hz.
%   high - high limit on passband, in Hz.
%   f   -   Sampling frequency (Hz)
%This function uses a BINARY window on the fourier transform.
%%Analysis without autocorrelation - just to find the modulation
Y = fft(A);
L = size(A);
x = f/2 * linspace(0, 1, L(2) / 2);
binwin = (x > low) & (x < high);
binwin = [binwin fliplr(binwin)];
Yw = binwin .* Y;
c = ifft(Yw); 
sq = square(unwrap(angle(hilbert(A))));
chop_on = (sq == 1);
chop_off = (sq == -1);
%TODO: Add in the phase information by using the MAXLAG feature of xcorr.
