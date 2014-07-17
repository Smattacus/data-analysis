function l = lorentz(x, x0, I, gamma)
%Quick function to calculate a given Lorentzian corresponding to the x data
%array. The function used is the following convention:
%
%L(x, x0 I, gamma) = I * gamma^2 / (gamma^2 + (x - x0)^2)
%
l = I * gamma^2 ./ (gamma^2 + (x - x0).^2);