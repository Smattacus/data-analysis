%%
%Script to average the autocorrelation of the GAGE data in the LIF
%acquisitions in a given directory
dyename = 'dnhere';
addpath('/Users/smattingly/Programs/MATLAB/');
cd('/Users/smattingly/Data/11-5-2013/V2PT/');
cd(strrep(dyename, '611_', '611.'));
%%
%Read in the files, take their autocorrelation, and average the spectra
list = dir('*.h5');
Nf = size(list,1);
temp = h5read(list(1).name, '/GAGE_CHAN_B');
N = size(temp,1);
acorrs = zeros(Nf, 2 * N - 1);
for i=1:Nf
    %Order doesn't matter, so we're just going to average everything.
    %Don't worry about contsructing filenames.
    fn = list(i).name;
    data = h5read(fn, '/GAGE_CHAN_B');
    ac = xcorr(double(data), 'unbiased');
    acorrs(i, :) = ac;
end
gagemeans = mean(acorrs);
%%
%Save the stuff:
save(sprintf('/Users/smattingly/Data/11-5-2013/V2PT/GageMeans/%s_gagemean.mat', dyename), 'gagemeans');
exit;