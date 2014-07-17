function [Vp, n, Te] = ReadParams(fileToRead1)
%
% [Vp, n, Te] = ReadParams(inputFile)
%
%ReadParams
%  Reads the density, plasma potential, and electron temperature.
%  FILETOREAD1:  file to read
% [n, Vp, Te]:   Output array of density, phi, and electron temp
% respective.

%  Auto-generated by MATLAB on 23-Feb-2013 17:38:30

DELIMITER = '=';
HEADERLINES = 2;

% Import the file
newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
Vp = newData1.data(1);
n = newData1.data(2);
Te = newData1.data(3);