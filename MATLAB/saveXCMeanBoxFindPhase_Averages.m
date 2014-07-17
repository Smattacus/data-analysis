function saveXCMeanBoxFindPhase_Averages(path, savename, sq_freq, acq_freq)
%This function generates and saves the averaged cross laser correlation
%function for a given folder.
%
%This functioni also saves the mean LIF phtoon count per file in a separate
%array.
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
noise1 = sqrt(mean(T1) + mean(B1));
noise2 = sqrt(mean(T2) + mean(B2));
ave1 = mean(D1);
ave2 = mean(D2);
N = size(D1,1);
xcs = zeros(size(D1, 1), size(D1,2) * 2 - 1);
xcs_lifscaled = xcs;
ac_ch1 = xcs;
ac_ch2 = xcs;
for i=1:size(D1, 1);
   D1(i,:) = D1(i,:) - mean(D1(i,:));
   D2(i,:) = D2(i,:) - mean(D2(i,:));
   xcs(i,:) = xcorr(D1(i,:), D2(i,:), 'unbiased');
   ac_ch1(i,:) = xcorr(D1(i,:), 'unbiased');
   ac_ch2(i,:) = xcorr(D2(i,:), 'unbiased');
   xcs_lifscaled = xcorr(D1(i,:) / ave1(i), D2(i,:) / ave2(i), 'unbiased');
end
xcmean = (mean(xcs));
xc_std = std(xcs);
xc_lifscaled_mean = mean(xcs_lifscaled);
xc_lifscaled_std = std(xcs_lifscaled);
ac_ch1_mn = mean(ac_ch1);
ac_ch2_mn = mean(ac_ch2);
ac_ch1_std = std(ac_ch1);
ac_ch2_std = std(ac_ch2);

save(savename, 'xcmean', 'xc_std', 'ac_ch1_mn', 'ac_ch2_mn', 'ac_ch1_std', ...
'ac_ch2_std', 'ave1', 'ave2', 'xc_lifscaled_mean', 'xc_lifscaled_std');