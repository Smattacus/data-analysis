function saveXCMeanBoxFindPhase_SumDiff(path, savename, sq_freq, acq_freq)
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
clear T1 B1 T2 B2
N = size(D1,1);
xcs = zeros(1, size(D1,2) * 2 - 1);
ac_ch1 = xcs;
ac_ch2 = xcs;
for i=1:size(D1, 1);
   D1(i,:) = D1(i,:) - mean(D1(i,:));
   D2(i,:) = D2(i,:) - mean(D2(i,:));
   xcs = xcs + xcorr(D1(i,:), D2(i,:), 'unbiased');
   ac_ch1 = xcs + xcorr(D1(i,:), 'unbiased');
   ac_ch2 = xcs + xcorr(D2(i,:), 'unbiased');
end
xcmean = xcs / size(D1,1);
xc_std = std(xcs);
ac_ch1_mn = ac_ch1 / size(D1,1);
ac_ch2_mn = ac_ch2 / size(D1,1);
ac_ch1_std = std(ac_ch1);
ac_ch2_std = std(ac_ch2);

save(savename, 'xcmean', 'xc_std', 'ac_ch1_mn', 'ac_ch2_mn', 'ac_ch1_std', 'ac_ch2_std');