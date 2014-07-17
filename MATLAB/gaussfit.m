function [A, mu, var] = gaussfit(x, y, h)
%Function for performing a gaussian fit of an array of data. Uses the
%builtin MATLAB function polyfit.
%
% [A, mu, var] = gaussfit(x, y)
%
%This function takes the log of the y values, so it is recommended that
%these values all be positive. This returns values according to
%
%f(x) = A exp(-(x-mu)^2 / (2 * var));
%
%
%INPUTS:
%   x   - Abscissa.
%   y   - Ordinate to be fitted.
%   h   - (0 - 1) Scaling parameter for point rejection. All points below this
%   fraction of the maximum will be rejected. A good start is 0.2.
%
%OUTPUTS:
%   A       - Amplitude of exponential term.
%   mu      - Center of gaussian.
%   var   - Variance

%Cut out the terms that are too low:
fh = find(y < (max(y) * h));
x(fh) = [];
y(fh) = [];
%Center the x - axis according to the mean of the terms (to be corrected by
%the fit)
m = mean(x);
p = polyfit(x - m, log(y), 2);
p
var = - 1 / (p(1) * 2);
mu = p(2) * 2 * var;
A = exp(p(3) + mu^2 / (2 * var));
mu = mu + m;
