function [f, p] = fft_plot_file(filename, dt);
%Quick function to plot a gagescope data file
f1 = fopen(filename);
A = fread(f1, 'short=>double');
A = A - mean(A);
[f p] = spec(A, dt);
semilogy(f, abs(p));
