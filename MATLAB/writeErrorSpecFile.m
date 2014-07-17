function writeErrorSpecFile(filename, tc)
%This is a function to write the spectra and error bars for said spectra
%using an exponential windowing function to a file.
%
%writeErrorSpecFile(filename, tc)
%
%INPUTS:
%   filename    - Filename of the matlab .mat file. Given like
%   "filename.mat"
%   tc          - Correlation time of the correlations.
%This function just reads in the .mat file denoted by filename, windows the
%cross correlations and standard deviations, and writes the spectra and
%error spectra to filename_err.mat
%
%This function assumes the correlations are named like:
%xcmean
%xc_std
%ac_ch1_mn
%ac_ch1_std
%ac_ch2_mn
%ac_ch2_std
load(sprintf('%s.mat', filename));
N = (size(xcmean,2) + 1) / 2;
t = (-(N-1):(N-1))/1e5;
win = exp(-(t/tc).^2/2);
[f, gxmn] = spec(xcmean .* win, 1e-5);
[f, gxsd] = spec(xc_std .* win, 1e-5);
[f, ga1mn] = spec(ac_ch1_mn .* win, 1e-5);
[f, ga2mn] = spec(ac_ch2_mn .* win, 1e-5);
[f, ga1sd] = spec(ac_ch1_std .* win, 1e-5);
[f, ga2sd] = spec(ac_ch2_std .* win, 1e-5);
save(sprintf('%s_err.mat', filename), 'f', 'gxmn', 'gxsd', 'ga1mn', 'ga2mn', 'ga1sd', 'ga2sd');