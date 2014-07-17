function [f, spec_matrix] = getSpecMatrixXcorr(directory, card, column1, column2, freq, freq_lim)
%Function to find the matrix of spectra in time corresponding to a set
%of control feedback readings taken over time.
%
%This is meant to be used with the controllooger.vi program, and reads the
%output files from that.
%
%INPUTS:
%   directory - Location of the files.
%   card      - Search parameter to find certain files in the directory.
%                   Can be left as ''. The search string is
%                   'Directory/*<card>*.txt'
%   column    - Which column of data to use.
%
%OUTPUTS:
%   spec_matrix - 2D matrix of spectra vs time.
%
%For now, the column of data for each file is autocorrelated and then
%windowed according to a correlation time of 0.01s.

tc = 0.01;
files = dir(sprintf('%s/*%s*.txt', directory, card));
f = fopen(files(1).name);
data = textscan(f, '%f %f %f %f %f %f', 'HeaderLines', 23);
fclose(f);
d = data{column1};
N = size(d,1);
t = (-(N-1):(N-1))/freq;
win = exp(-(t/tc).^2/2);
if freq_lim > 0
    [ft, gt] = spec(xcorr(d, 'unbiased'), 1/freq);
    ig = find(abs(ft) < freq_lim);
    spec_matrix = zeros(size(files,1), size(ig,2));
else
    spec_matrix = zeros(size(files,1), size(d,1) * 2 - 1);
end
for i=1:size(files,1)
    temp = dir(sprintf('*num%i.txt', i-1));
    f = fopen(temp.name);
    data = textscan(f, '%f %f %f %f %f %f', 'HeaderLines', 23);
    fclose(f);
    d1 = data{column1};
    d1 = d1 - mean(d1);
    d2 = data{column2};
    d2 = d2 - mean(d2);
    ac = xcorr(d1, d2, 'unbiased');
    [f, gw] = spec(ac .* win.', 1/freq);
    if freq_lim > 0 
        spec_matrix(i,:) = gw(ig);
    else
        spec_matrix(i,:) = gw;
    end
end
if freq_lim > 0 
    f = f(ig);
end