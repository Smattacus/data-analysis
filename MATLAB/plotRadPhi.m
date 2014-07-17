function A = plotRadPhi(path, suffix, offset)
%Generates a plot of density vs angle of a plasma. This takes into account
%the length of the langmuir probe, and ASSUMES that the probe passes
%through the center of the plasma. See below for details.
%
% plotRadPhi(path, suffix);
%
%INPUTS:
% path      - Directory holding the file of interest.
% suffix    - The UNIQUE suffix to the filename.
%
%OUTPUTS:
% DATA      - Matrix of [radius, Te, n, Vp] in Nx4
%
%Note about suffixes: For the filename
%
%Feb-20-2013-01-05-48-PM-0,0--25,200-full_22W_30_param.txt
%
%the suffix is "full_22W_30_param.txt". Each data run should have its own
%unique suffix, and this function will plot all the points from a given
%suffix.
%
%
%In the case the probe passes through the center of the plasma, 
%it traces out an
%arc which can be used to find x and y and thus (r, theta) in a polar
%coordinate system centered at the center of the plasma. In this case, 
%x = R(1-cos(theta)) and y = R sin(theta), where theta=0 is the probe tip
%at the center of the plasma and R is the length of the probe.
%x and y are the coordinate in the new coordinate system so 
%r = sqrt(x^2 + y^2)

R = 10.5; %centimeters based off Peter Haugen's measurement.
A = getAngDensity(path, suffix);
A(1,:) = A(1,:) + offset;
x = R* (1 - cos(pi / 180 * A(1,:)));
y = R * sin(pi / 180 * A(1,:));
rvec = sqrt(x.^2 + y.^2);
angles = A(1,:);
rvec = rvec .* (angles > 0) - rvec .* (angles < 0);
plot(rvec, A(2,:));
xlabel('Radius (cm)');
ylabel('Phi (volts)');
suffix = strrep(suffix, '_', '-');
tn = sprintf('Radius vs phi plot for LMP Run %s', suffix);
title(tn);