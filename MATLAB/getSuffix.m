function suffix = getSuffix(filename)
%
% angle = getAngle(filename)
%
%Helper function to get the suffix from a langmuir probe filename which has
%been named with Feng's naming convention.
%Inputs:
%   filename - Input filename. Will work with time series or LMP data.
%Outputs:
%   suffix    - suffix of the file, minus the .txt.
[s, r] = strtok(filename, ',');
[s, r] = strtok(r, ',');
[s, r] = strtok(r, '-');
suffix = strrep(r, '-', '');
end
    


