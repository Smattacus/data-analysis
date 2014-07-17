function [norm_str_1, norm_str_2] = get_LIFNormStrength(A, left, right, tc)
%Function to give the normalized, integrated LIF strength
%for a given LIF file.
%The routine reads a file, sums and subtracts off the mean, autocorrelates
%the two arrays, then performs a windowed FFT according to tc.
%
%[o1, o2] = get_LIFNormStrength(A, left, right, tc)
%
%INPUTS:
%   A       - Input 32xN array.
%   left    - Left bound to integrate over in f(v) space.
%   right   - Right bound.
%   tc      - Correlation length to window by.
%
%OUTPUTS:
%   o1      - Output from channels 1 - 16.
%   o2      - Output from channels 17 - 32.
%
s1 = sum(A(1:16, :));
m1 = mean(s1);
s1 = s1 - m1;
s2 = sum(A(17:32,:));
m2 = mean(s2);
s2 = s2 - m2;
N = size(s1, 2);
t = (-(N-1):(N-1))/1000000;
%x1 = xcorr(s1, 'unbiased');
%x2 = xcorr(s2, 'unbiased');
[f, g1] = spec(s1, 1e-6);
[f, g2] = spec(s2, 1e-6);
nmid = (size(g1,2) + 1)/2;
norm_str_1 = 2 * max(abs(g1(left:right))) / m1;
norm_str_2 = 2 * max(abs(g2(left:right))) / m2;
end