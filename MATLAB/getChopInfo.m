function [fns centers deltas] = getChopInfo(filename)
%
%[fn centers deltas] = getChopInfo(filename)
%
%Loads the information from a chop info file of the form:
%
%FILENAME   CHOP_FREQUENCY_CENTER   WINDOW_WIDTH
%<val>      <value>                 <value>
%
%This simply reads out the above file and returns the columns of filenames,
%chop freq centers, and widths for use in genPhase().
%
format = '%s';
f = fopen(filename, 'r');
strs = textscan(f, format, 'Delimiter', '', 'EmptyValue', NaN, 'HeaderLines', 1, 'ReturnOnError', false);
strs = strs{:};
centers = zeros(size(strs,1), 1);
deltas = zeros(size(strs,1), 1);
fns = cell(size(strs,1), 1);
for i=1:size(strs,1)
    [temp, r] = strtok(strs{i});
    fns{i} = temp;
    [temp, r] = strtok(r);
    centers(i) = str2num(temp);
    [temp, r] = strtok(r);
    deltas(i) = str2num(temp);
end
end