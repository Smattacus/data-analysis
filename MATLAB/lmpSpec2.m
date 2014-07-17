function [freq, Q, H, data] = lmpSpec2(filename, filterWidth, F)
%[freq, Q, H, data] = lmpSpec2(filename, filterWidth, F)
%
%Generate a windowed 2D Fourier transform of a langmuir probe data set.
%
%Inputs:
%filename       - String of the filename of the LMP array file.
%filterWidth    - Width of the Gaussian window. Should be longer than
%correlation time. For now, code coerces to .015. (line 42 abouts)
%F              - Sampling frequency for the probes.
%Outputs:
%freq             - Frequency array
%Q              - Log normalized fftshifted power spectrum
%H              - Raw fft2 spectrum
%data           - Data cell array read from the file. Constructed as
%                 {LP; AP1; AP2; ...; AP8}


data = getLMPArray(filename);
whos data
AP1 = data{:,2};
AP2 = data{:,3};
AP3 = data{:,4};
AP4 = data{:,5};
AP5 = data{:,6};
AP6 = data{:,7};
AP7 = data{:,8};
AP8 = data{:,9};
G = [AP1 AP2 AP3 AP4 AP5 AP6 AP7 AP8];
whos G
s = size(G);
clearvars AP1 AP2 AP3 AP4 AP5 AP6 AP7 AP8
nyq = F / 2;
N = s(1) / 2;
freq = nyq * ((-N):(N-1)) / N;
mG = mean(G);
G = G - ones(s(1),1) * mG; %Subtract off mean
sG = std(G); %RMS fluctuation
%Later we will try some SVD bsnznsnznsnz (banana business)
%For now, do some xcorrs between the probes
c1 = xcorr(G(:,1), 'unbiased');
c2 = xcorr(G(:,1), G(:,2), 'unbiased');
c3 = xcorr(G(:,1), G(:,3), 'unbiased');
c4 = xcorr(G(:,1), G(:,4), 'unbiased');
c5 = xcorr(G(:,1), G(:,5), 'unbiased');
c6 = xcorr(G(:,1), G(:,6), 'unbiased');
c7 = xcorr(G(:,1), G(:,7), 'unbiased');
c8 = xcorr(G(:,1), G(:,8), 'unbiased');
t = (-(s(1) - 1):s(1) - 1) / F; %Create the time axis.
tc = filterWidth; %Filter width (Longer than correlation length 
                   %of the phenomena of interest)
%Make a windowed cross-corr array of the probe array.
C = [exp(-(t'/tc).^2/2).*c1, exp(-(t'/tc).^2/2).*c2, ...
    exp(-(t'/tc).^2/2).*c3, exp(-(t'/tc).^2/2).*c4, ...
    exp(-(t'/tc).^2/2).*c5, exp(-(t'/tc).^2/2).*c6, ...
    exp(-(t'/tc).^2/2).*c7, exp(-(t'/tc).^2/2).*c8];
H=fft2(C(N + (1:s(1)),:) - ones(s(1),1)*mean(C(N + (1:s(1)),:)));
Q=fftshift(log(abs(H).^2));
end