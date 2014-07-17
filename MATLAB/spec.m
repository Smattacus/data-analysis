function [f ,P]=spec(g,dt);
%Discrete Fourier Transform that conserves signal energy.
%Result is shifted so that f=0 is in the center
%  [f,P]=spec(g,dt)  Input function g(t) and the sample spacing dt
%
%  Returns the transform P (complex) and the array of frequencies f
n=length(g);
f0=1/abs(n*dt);
P=fft(g)/sqrt(n);
if rem(n,2)==0
    f=f0*[0:((n-2)/2) (-(n/2)):-1];
else
    f=f0*[0:((n-1)/2) (-(n-1)/2):-1];
end
f=fftshift(f);
P=fftshift(P);
return