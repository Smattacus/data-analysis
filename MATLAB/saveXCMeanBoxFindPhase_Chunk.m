function saveXCMeanBoxFindPhase_Chunk(path, savename, sq_freq, acq_freq)
%This function generates and saves the averaged cross laser correlation
%function for a given folder.
%
% saveXCMean(path, savename, phases)
%
% INPUTS:
% path          - Directory to save the *.mat file to.
% savename      - Will be used to create the .mat file, according to
%                   savename.mat
% phases        - [PMT1, PMT2] array of phases to use.
%
cd(path);
%[T1, B1, T2, B2] = genAllTB(path, info_file, parallel, nworkers);
%For now, use defaults:
[T1, B1, T2, B2] = genAllTPBoxFindPhase(path, sq_freq, acq_freq);
D1 = T1 - B1;
D2 = T2 - B2;
N = size(D1,1);
xcs = zeros(size(D1, 1), size(D1,2) * 2 - 1);
ac_ch1 = xcs;
ac_ch2 = xcs;
for i=1:size(D1, 1);
   D1(i,:) = D1(i,:) - mean(D1(i,:));
   D2(i,:) = D2(i,:) - mean(D2(i,:));
   xcs(i,:) = xcorr(D1(i,:), D2(i,:), 'unbiased');
   ac_ch1(i,:) = xcorr(D1(i,:), 'unbiased');
   ac_ch2(i,:) = xcorr(D2(i,:), 'unbiased');
end
xcmean = (mean(xcs));
xc_std = std(xcs);
ac_ch1_mn = mean(ac_ch1);
ac_ch2_mn = mean(ac_ch2);
ac_ch1_std = std(ac_ch1);
ac_ch2_std = std(ac_ch2);
xcmean_50 = mean(xcs(1:50,:));
xcmean_60 = mean(xcs(1:60,:));
xcmean_70 = mean(xcs(1:70,:));
xcmean_80 = mean(xcs(1:80,:));
xcmean_90 = mean(xcs(1:90,:));
xcmean_100 = mean(xcs(1:100,:));
xcmean_110 = mean(xcs(1:110,:));
save(savename, 'xcmean', 'xc_std', 'ac_ch1_mn', 'ac_ch2_mn', 'ac_ch1_std', 'ac_ch2_std', 'xcmean_50', 'xcmean_60', 'xcmean_70', 'xcmean_80', 'xcmean_90', 'xcmean_100', 'xcmean_110');