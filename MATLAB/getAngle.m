function angle = getAngle(filename)
%
% angle = getAngle(filename)
%
%Helper function to get the angle from a langmuir probe filename which has
%been named with Feng's naming convention.
%Inputs:
%   filename - Input filename. Will work with time series or LMP data.
%Outputs:
%   angle    - Angle in decimal degrees.
[s, r] = strtok(filename, ',');
t = strfind(r, '-');
%if TRUE then angle measure is positve, else negative.
if r(t(1)) == r(t(1) + 1)
   [t, r] = strtok(r, '-');
   [t, r] = strtok(r, '-');
   angle = - str2double(strrep(t, ',', '.'));
else
   [t, r] = strtok(r, '-'); 
   [t, r] = strtok(r, '-');
   angle = str2double(strrep(t, ',', '.'));
end
    


