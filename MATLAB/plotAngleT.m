function A = plotAngleT(path, suffix, s_b)
%Generates a plot of density vs angle of a plasma. This has not yet been
%normalized for the langmuir probe's length.
%
% plotAngleT(filename, suffix, scatter);
%
%INPUTS:
% path      - Directory holding the file of interest.
% suffix    - The UNIQUE suffix to the filename.
%scatter    - Boolean. true makes a scatter plot.
%
%OUTPUTS:
% DATA      - Matrix of [angle, Te, n, Vp] in Nx4
%
%Note about suffixes: For the filename
%
%Feb-20-2013-01-05-48-PM-0,0--25,200-full_22W_30_param.txt
%
%the suffix is "full_22W_30_param.txt". Each data run should have its own
%unique suffix, and this function will plot all the points from a given
%suffix.

A = getAngDensity(path, suffix);
if s_b
    scatter(A(1,:), A(2,:),'x');
else
    plot(A(1,:), A(2,:));
end
xlabel('Angle');
ylabel('Temperature (eV)');
suffix = strrep(suffix, '_', '-');
tn = sprintf('Angle vs Temperature plot for LMP Run %s', suffix);
title(tn);