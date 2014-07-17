function [f, p] = analyze_gage(filename, f_s, f_c)
%Function to analyze a gagescope data file.
%Takes the hilbert transform, subtracts the linear offset,
%and finds the spectrum of the remaining phase information.
%Inputs:
%   filename - filename of gagecard file.
%   f_s         - Sampling frequency
%   f_c         - Carrier frequency
f = fopen(filename);
x = fread(f, 'short=>double');
%Subtract offset
xn = x - mean(x);
ht = hilbert(xn);
l = size(ht);
l = l(1);
upha = unwrap(angle(ht));
%upha is the unwrapped phase angle:
%W * t + g + sin(v * t + p)
%we want the info in the sin, so remove W*t by using the
%carrier frequency information:
os = l * f_c / f_s * (1:1 / f_s:l / f_s);
uphas = upha - os; %Subtract linear offset
