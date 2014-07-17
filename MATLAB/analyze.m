%Test script for analyzing a gagescope file,
%with 4kHz chopping on the signal.
f = fopen('Tue_May_01_00h42m38s2012__plasma_LIF_1MHz_0_A.DAT');
A = fread(f, 'short=>double');
Fs = 1e6; %1 MHz sampling frequency
raw = A.';
A = A - mean(A);
xca = xcorr(A);
s = size(xca);
win = gausswin(s(1), 10);
wxc = win .* xca;
[fxcw, Pxcw] = spec(wxc, 1/f);
%Do some binary filtering - just take the peak close to 4khz
%Create a series of 1's near -4khz and 4 khz (+/- 25 Hz)
t = 0:2.097152:1e-6;

%%Analysis without autocorrelation - just to find the modulation
Y = fft(A);
L = max(size(A));
x = Fs/2 * linspace(0, 1, L / 2);
t = 1e-6 * L * linspace(0, 1, L); 
