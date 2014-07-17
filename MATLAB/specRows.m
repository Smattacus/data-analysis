function [f ,P]=specRows(g,dt)
%Discrete Fourier Transform that conserves signal energy.
%Result is shifted so that f=0 is in the center
%  [f,P]=specRows(g,dt)  Input function g(t) and the sample spacing dt
%
%  Returns the transform P (complex) and the array of frequencies f
%
% INPUTS:
%g      - NxM matrix of data. FFT will be performed across the N
%dimensions.
n=length(g);
f0=1/abs(n*dt);
P=fft(g,[], 1)/sqrt(n);
if rem(n,2)==0
    f=f0*[0:((n-2)/2) (-(n/2)):-1];
else
    f=f0*[0:((n-1)/2) (-(n-1)/2):-1];
end
f=fftshift(f);
P=fftshift(P,1);
return