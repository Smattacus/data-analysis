function [t,g]=ispec(P,df)
%ispec is the inverse of the Fourier Transform given by spec.m which is my 
%version that conserves signal energy and returns a frequency axis.
%Given the FFT P and the increment in frequency df:
%
%[t,g]=ispec(P,df) 
%
%   returns the time domain signal g and the time axis t.
%   In the place of df you can put the array f given by spec
if length(df)>1
    df=df(2)-df(1);
end
T=1/abs(df);
n=length(P);
g=ifft(ifftshift(P))*sqrt(n);
t=T*(1:n)/n;
return