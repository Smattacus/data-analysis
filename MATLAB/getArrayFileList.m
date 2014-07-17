function list = getArrayFileList(varargin)
%Function to get a list of the files corresponding to a particular data
%acquisition using the LMP probe array. See below for detailed
%instructions.
%
%list = getArrayFileList(path, suffix)
%list = getArrayFileList()
%
%Inputs:
% path      - Path to the directory of interest
% suffix    - Unique suffix for a given data run. See below.
% (none)    - Opens a dialog asking to choose a file.
%Outputs:
% list      - A cell array of filenames that fit the criteria.
%
%The filenames MUST be of the format
%
%<DATE>-<POSITION>-<ANGLE>-<SUFFIX>.TXT
%e.g.
%Feb-20-2013-01-05-48-PM-0,0--25,200-<suffix>.txt
%
%Upon giving a path and suffix, this routine will return all the files
%in the directory of path with the given suffix.
%

if nargin ~= 2
    disp('Select any file from a langmuir probe angular acquisition run');
    [fn, path, index] = uigetfile('Select a valid LMP param .txt file.', '*.txt');
    o = sprintf('%s has been selected', fn); disp(o);
    %TODO: This can be automated by just parsing by comma, then selecting
    %the second element and parsing by hyphen. Everything from last hyphen
    %to end is what's needed.
    disp('The suffix is after the "-" after the angle measurement to the end');
    disp('Include the file extension');
    r = input('Input unique suffix to this file: (e.g. etc_param.txt)', 's');
else
    path = varargin{1};
    r = varargin{2};
end

query = sprintf('%s/*%s', path, r);
list = dir(query);
end

