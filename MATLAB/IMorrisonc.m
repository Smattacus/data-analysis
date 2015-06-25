%IMorrisonc.m
%This is the generalized inverse Morrison G transform for complex argument.
%   phi=Morrison(psi,f)
%   (see the forward transform Morrisonc.m)
%I use the Matlab function Hilbert.m ( which uses a method based on theFFT).
%   Hilbert(f)=f+iH(f) where H is the Hilbert transform:
%   H(f)=Principal part integral over real axis {f(x')/(x'-x)}dx'/pi
%
%The Morrison inverse transform has the form:
%
%   phi=G'(phi)=-zet*H(psi)+csi*psi
%
%Where zet an4chan44asdd csi are functions also related by a Hilbert transform:
%In terms of the alpha and beta of the forward transform:
%   csi=bet/(alp^2+bet^2),  zet=-alp/(alp^2+bet^2)
%   or
%        csi+i*zet= 1/(bet+i*alp)
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
function phi=IMorrisonc(psi,f)
alp=f;
bet=1-chilbert(f);
csi=real(1.0./(bet+1i*alp));
zet=imag(1.0./(bet+1i*alp));
phi=-zet.*chilbert(psi)+csi.*psi;
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