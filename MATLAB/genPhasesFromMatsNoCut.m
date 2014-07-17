function [f, sxw, phases, ph_errs, mfile_list] = genPhasesFromMatsNoCut(directory, tc)
%This is a general function for determining the spectra, phases, and the 
%errors on the phases from all the *.mat files of a given directory. 
%
%[sxw, phases, ph_errs, mfile_list] = genPhasesFromMats(directory, Nc, tc,)
%
%
%Inputs:
% directory     - Path to the directory of .mat files.
% tc            - Width of windowing function, in seconds.
%
%Outputs:
% sxw           - Array of the spectra from the files. Goes like [Nf, N]
%                   where NF is the number of files.
% phases        - Phases from the spectra. Organized the same way as sxw.
% ph_errs       - Errors of the phases. Organized the same way as sxw.
% mfile_list    - List of mat files used in the spectra. Nf long,
%                  corresponds to the Nf dimension of the sxw array.
cd(directory);
mfile_list = dir (sprintf('%s/*.mat', directory));
load(mfile_list(1).name);
N = (size(xcsmean_lifscaled, 2) + 1) / 2;
Nf = size(mfile_list, 1);
t = (-(N-1):(N-1))/1e5;
win = exp(-(t/tc).^2/2);
sxw = zeros(Nf, 2 * N - 1);
phases = zeros(Nf, 2 * N - 1);
ph_errs = zeros(Nf, 2 * N - 1);
for i=1:size(mfile_list)
    load(mfile_list(i).name);
    temp = xcsmean_lifscaled - mean(xcsmean_lifscaled);
    ec = mean(std(xc_std_lifscaled));
    e_redc = ec * sqrt(sum(win.^2) / (2 * (2 * N - 1))) / sqrt(size(ave1,2));
    [f, sxw(i, :)] = spec(temp .* win, 1e-5);
    ph_errs(i,:) = e_redc ./ abs(sxw(i,:));
end
phases = angle(sxw);