function saveXCMean(path, info_file, savename, parallel, nworkers)
%This function generates and saves the averaged cross laser correlation
%function for a given folder.
%
% saveXCMean(path, info_file, savename, parallel, nworkers)
%
% INPUTS:
% path          - Directory to save the *.mat file to.
% info_file     - File containing the filenames, chopinfo, 
% savename      - Will be used to create the .mat file, according to
%                   savename.mat
% parallel      - Boolean to use parallelization.
% nworkers      - Number of client workers to open.
%
cd(path);
if size(info_file,1) == 0;
    info_file = dir('*chopinfo.txt');
    info_file = info_file(1).name;
end
%[T1, B1, T2, B2] = genAllTB(path, info_file, parallel, nworkers);
%For now, use defaults:
[T1, B1, T2, B2] = genAllTB(path, info_file);
D1 = T1 - B1;
D2 = T2 - B2;
N = size(D1,1);
xcs = zeros(size(D1, 1), size(D1,2) * 2 - 1);
ac_ch1 = zeros(size(D1, 1), size(D1,2) * 2 - 1);
ac_ch2 = zeros(size(D1, 1), size(D1,2) * 2 - 1);
for i=1:size(D1, 1);
   D1(i,:) = D1(i,:) - mean(D1(i,:));
   D2(i,:) = D2(i,:) - mean(D2(i,:));
   xcs(i,:) = xcorr(D1(i,:), D2(i,:), 'unbiased');
   ac_ch1(i,:) = xcorr(D1(i,:), 'unbiased');
   ac_ch2(i,:) = xcorr(D2(i,:), 'unbiased');
end
xcmean = (mean(xcs));
xcstd = std(xcs);
ac_ch1_mn = mean(ac_ch1);
ac_ch1_std = std(ac_ch1);
ac_ch2_mn = mean(ac_ch2);
ac_ch2_std = std(ac_ch2);
save(savename, 'xcmean', 'xcstd', 'ac_ch1_mn', 'ac_ch1_std', 'ac_ch2_mn', 'ac_ch2_std');