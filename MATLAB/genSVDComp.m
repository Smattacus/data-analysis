function [angle, Gsin, Gcos] = genSVDComp(filename)
%Function to generate an SVD G vector where G is the ansantz
%G ~ x(t)f(theta). This is to be used in a larger array
%which is used in SVD.
%
%Inputs:
% filename      - Desired array probe data filename
%Outputs:
% angle         - Angle corresponding to the filename.
% G             - G vector function

data = getLMPArray(filename);
angle = getAngle(filename);
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
n2 = size(x, 2);
%Take the fft along the 8 probes (but not in time)
spatmod = fft(x.').';
%Choose power conserving normalization
spatmod = spatmod / sqrt(n2);
%Select the m=1,-1 mode.
m1_t = spatmod(:,6);
m1_nt = spatmod(:,4);
%Add/diff m=1 and m=-1 parts for cosine/sine terms. 
mode_cos = (m1_t + m1_nt)/2;
mode_sin = (m1_t - m1_nt) / (2*1i);
mcos_mean = mean(mode_cos);
msin_mean = mean(mode_sin);
Gcos = xcorr(mode_cos - mcos_mean, LP - mean(LP), 'unbiased');
Gsin = xcorr(mode_sin - msin_mean, LP - mean(LP), 'unbiased');
end