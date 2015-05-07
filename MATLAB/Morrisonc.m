%Morrisonc.m
%This is the generalized Morrison G transform for complex argument.
%This version also works for real argument but requires the function chilbert.
%   Psi=Morrison(phi,f)
%   (see the inverse transform IMorrisonc.m)
%I use the Matlab function Hilbert.m ( which uses a method based on theFFT).
%   Hilbert(f)=f+iH(f) where H is the Hilbert transform:
%
%   H(f)=principal part integral over real axis {f(x')/(x'-x)}dx'/pi
%
%The Morrison transform has the form:
%
%   Psi=G(phi)=-alp*H(phi)+bet*phi
%
%Where alp and bet are functions also related by a Hilbert transform:
%
%   bet=-H(alp)+bet0  (bet0 is actually the value of beta at argument =inf)
%
%We normalize the transform so that
%   1=-H(alp)+bet     (in otherwords, we make bet0=1)
%The transform breaks down if abs(i*alp+bet)=0 somewhere on the real axis.
%This normalization prevents the zero from happening.
%
%We use the real input function f to produce alp and bet
%  alp=f
%  bet=-H(f)+1
%
%The motivation is the expansion in Case-Van Kampen modes
%In this case, for the ion acoustic wave,
%  f=-pi*Cs^2*df/fv=Imag(epsilon)
%In the case of the electron plasma wave
%  f=-pi*(wp/k)^2*df/dv
%In both cases bet=Real(epsilon)
%
%Old routines called MGT and IMGT
%This version written 9/2006
%************************************************************
%
function psi=Morrisonc(phi,f)
alp=f;
bet=1-chilbert(f);
psi=-alp.*chilbert(phi)+bet.*phi;
return

%chilbert.m
%hilbert transform of complex argument
%note that for real argument the result is not the same as 
%using Matlab hilbert.
function z=chilbert(f)
if isreal(f)
    z=imag(hilbert(f));
else
    z=imag(hilbert(real(f)))+1i*imag(hilbert(imag(f)));
end
return