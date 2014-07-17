function [A] = getAngDensity(varargin) 
%This function returns a matrix which gives the [angle, temp, density, phi]
%for a given set of files.
%
%USAGE:
%NO INPUTS: 
%The function will prompt the user for an example filename
%and a suffix particular to the data for that the file. For example, the filename
%Feb-20-2013-01-05-48-PM-0,0--25,200-full_22W_30_param.txt
%has the suffix "full_22W_30_param.txt" which can be used to find all the
%other points which were taken in this data run.
%
%MANUAL INPUTS:
%A = get_ang_density(path, suffix)
%Path       - Path to the folder holding the desired data run.
%Suffix     - Identifying Suffix for a given data run.
%
%
%In order to use this program, simply give it one of the files from the
%data run and it will attempt to find all the other files from that
%particular data run and generate the matrix of [angle, temp, density, phi]
%which is an Nx4 matrix.

if nargin ~= 2
    disp('Select any file from a langmuir probe angular acquisition run');
    [fn, path, index] = uigetfile('Select a valid LMP param .txt file.', '*.txt');
    o = sprintf('%s has been selected', fn); disp(o);
    %TODO: This can be automated by just parsing by comma, then selecting
    %the second element and parsing by hyphen. Everything from last hyphen
    %to end is what's needed.
    r = input('Input unique suffix to this file: (e.g. etc_param.txt)', 's');
else
    path = varargin{1};
    r = varargin{2};
end

query = sprintf('%s/*%s', path, r);
list = dir(query);

A = zeros(4, size(list, 1))
for i=1:size(list,1)
    cfile = list(1).name;
    cangle = getAngle(cfile);
    [cn, cVp, cTe] = ReadParams(cfile);
    A(:, i) = [cangle; cn; cVp; cTe];
end


end

