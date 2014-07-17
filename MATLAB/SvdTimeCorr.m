function [T, XC] = SvdTimeCorr(filename)
%%%Function which takes a given langmuir probe file array (filename),
%generates the SVD principal axes, and correlates these
%principal axes with the standalone probe.
%
%Inputs:
% filename      - Name of the langmuir probe array file
%
%Outputs:
% T        - Principal axes of the SVD. Nx8 array where N is the number of
%               samples.
% XC       - Cross correlation of principal axes with the standalone probe.
%               2Nx8.
%
[X V S U] = genSVDFromArray(filename);
data = getLMPArray(filename);
T = V * S;
lmp_alone = data{1};
XC = []; 
for i=1:8
    XC = [XC, xcorr(lmp_alone, T(:, i), 'unbiased')]; 
end