%First iteration script to try and find the SVD decomposition of the 
%m=1 mode. Other modes can be added easily.

[filename, path, fn] = uigetfile('*.txt');
filename = [path filename];

data = getLMPArray(filename);

LP = data{:,1};
AP1 = data{:,2};
AP2 = data{:,3};
AP3 = data{:,4};
AP4 = data{:,5};
AP5 = data{:,6};
AP6 = data{:,7};
AP7 = data{:,8};
AP8 = data{:,9};

F = 55.5e3;
%Phase dithering is unforgiving of fft
x = [AP1 AP2 AP3 AP4 AP5 AP6 AP7 AP8];
n1 = size(x, 1);
n2 = size(x, 2);
%2D FFT. Compensate for normalization in the mode direction.
%Add m=1 and m=-1 parts. 
%(The time direction will be removed with ifft)
spatmod = fft(x.').';
%Select the m=1 mode.
m1_t = spatmod(:,6) / sqrt(n2);
m1_mean = mean(m1_t);
G = xcorr(m1_t, LP);
