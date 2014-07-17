function qplot_FFT(A, tc, min_f, max_f)
%Quick function to plot the FFT with a typical window.
%Generates a plot after windowing the autocorrelation with an exponential
%with a correlation time of tc.
%
%qplot_FFT(A, tc, min_f, max_f)
%
%INPUTS:
%  A            - 32xN array of data.
%  tc           - Correlation time. Typically 0.01
%  max_f        - Maximum frequency to plot.
s1 = sum(A(1:16, :));
s2 = sum(A(17:32,:));
s1 = s1 - mean(s1);
s2 = s2 - mean(s2);
x1 = xcorr(s1, 'unbiased');
x2 = xcorr(s2, 'unbiased');
N = size(s1, 2);
t = (-(N-1):(N-1)) /1000000;
[f1, g1] = spec(x1 .* exp(-(t/tc).^2/2), 1e-6);
[f2, g2] = spec(x2 .* exp(-(t/tc).^2/2), 1e-6);
ig = find((f1 < max_f) - (f1 < min_f));
figure(1); semilogy(f1(ig), abs(g1(ig)));
xlabel('Frequency (Hz'); ylabel('Abs |P|'); title('Channels 1 - 16');
figure(2); semilogy(f2(ig), abs(g2(ig)));
xlabel('Frequency (Hz'); ylabel('Abs |P|'); title('Channels 17 - 32');
end