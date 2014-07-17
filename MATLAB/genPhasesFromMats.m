function [f, sxw, phases, ph_errs, mfile_list] = genPhasesFromMats(directory, Nc, tc, shift, mat_string)
%This is a general function for determining the spectra, phases, and the 
%errors on the phases from all the *.mat files of a given directory. 
%
%[sxw, phases, ph_errs, mfile_list] = genPhasesFromMats(directory, Nc, tc, shift)
%
%
%Inputs:
% directory     - Path to the directory of .mat files.
% Nc            - Number of points in the middle to use for a spectrum.
% tc            - Width of windowing function, in seconds.
% shift         - Number of elements to circular shift by.
% mat_string    - String to search for the mat files by. Queried with  
%                   dir('<mat_string>*.mat');
%
%Outputs:
% sxw           - Array of the spectra from the files. Goes like [Nf, N]
%                   where NF is the number of files.
% phases        - Phases from the spectra. Organized the same way as sxw.
% ph_errs       - Errors of the phases. Organized the same way as sxw.
% mfile_list    - List of mat files used in the spectra. Nf long,
%                  corresponds to the Nf dimension of the sxw array.
cd(directory);
mfile_list = dir (sprintf('%s/%s*.mat', directory, mat_string));
load(mfile_list(1).name);
N = (size(xcsmean_lifscaled, 2) + 1) / 2;
Nf = size(mfile_list, 1);
t = (-(N-1):(N-1))/1e5;
win = exp(-(t/tc).^2/2);
ui = int32(N - Nc/2);
hi = int32(ui + Nc);
sxw = zeros(Nf, size(ui:hi,2));
phases = zeros(Nf, size(ui:hi,2));
ph_errs = zeros(Nf, size(ui:hi,2));
for i=1:size(mfile_list)
    load(mfile_list(i).name);
    temp = xcsmean_lifscaled - mean(xcsmean_lifscaled);
    ec = mean(std(xc_std_lifscaled));
    e_redc = ec * sqrt(sum(win(ui:hi).^2) / (2 * Nc)) / sqrt(size(ave1,2));
    [f, sxw(i, :)] = spec(circshift(temp(ui:hi) .* win(ui:hi), [0,shift]), 1e-5);
    ph_errs(i,:) = e_redc ./ abs(sxw(i,:));
end
phases = angle(sxw);