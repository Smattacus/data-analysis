function data = getLMPArray(filename)
%Exactly as the name implies. Just reads out an LMP array data file.
%
%data = getLMPArray(filename)
%
%Inputs:
%   filename    - String of the desired filename to read.
%Outputs:
%   data        - data = [LP; AP1; AP2; AP3; AP4; AP5; AP6; AP7; AP8]
%The names of the elements in data are the same as that in a typical
%array data file.
%
% Import the file
% Import the file
format = '%f%f%f%f%f%f%f%f%f';
f = fopen(filename, 'r');
data = textscan(f, format, 'Delimiter', '\t', 'EmptyValue', NaN, ...
    'HeaderLines', 1, 'ReturnOnError', false);
end