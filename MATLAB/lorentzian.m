function f = lorentzian(b, x)
%Using the standard Physics three parameter Lorentzian distribution:
%
% f(x; x_0, gamma, I) = I / (1 + ((x - x_0) / gamma)^2)
%
% f = lorentzian(b, x)
%
%This function is written to satisfy the @modelfun type in nlinfit.
%The coefficients vector b corresponds to I, gamma, and x_0:
%
% b(1) = x_0
% b(2) = gamma
% b(3) = I
%
x_0 = b(1);
gamma = b(2);
I = b(3);
f = I ./ (1 + ((x - x_0)/gamma).^2);