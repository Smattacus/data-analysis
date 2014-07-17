function saveXCMeanBoxFindPhase_ErrProp(path, savename, sq_freq, acq_freq)
%This function generates and saves the averaged cross laser correlation
%function for a given folder.
%
%This function also saves the mean LIF photon count per file in a separate
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
ave1 = mean(D1.');
ave2 = mean(D2.');
N = size(D1,1);
xcs = zeros(size(D1, 1), size(D1,2) * 2 - 1);
xcs_lifscaled = xcs;
ac_ch1 = xcs;
ac_ch2 = xcs;
xcs_properr = xcs;
xcs_properr_scaled = xcs;
for i=1:size(D1, 1);
   display(sprintf('On run # %d in saveXCMeanBoxFindPhase_ErrProp', i));
   tic
   D1(i,:) = D1(i,:) - mean(D1(i,:));
   D2(i,:) = D2(i,:) - mean(D2(i,:));
   xcs(i,:) = xcorr(D1(i,:), D2(i,:), 'unbiased');
   xcs_lifscaled(i,:) = xcorr(D1(i,:) / ave1(i), D2(i,:) / ave2(i), 'unbiased');
   ac_ch1(i,:) = xcorr(D1(i,:), 'unbiased');
   ac_ch2(i,:) = xcorr(D2(i,:), 'unbiased');
   toc
   display(sprintf('Correlations Done, Run #%d', i));
   tic;
   xcs_properr(i,:) = xcorr_err(D1(i,:), D2(i,:), noise1(i), noise2(i));
   toc
   %Standard error prop for putting the rescaled data into the xcorr_err
   %function
   display(sprintf('Err 1 Done, Run #%d', i));
   tic
   xcs_properr_scaled(i,:) = xcorr_err(D1(i,:) / ave1(i), D2(i,:)/ave2(2), ...
       noise1(i) / ave1(i), noise2(i)/ave2(i));
   display(sprintf('Err 2 Done, Run #%d', i));
   toc
end
xcmean = (mean(xcs));
xcsmean_lifscaled = mean(xcs_lifscaled);
xc_std = std(xcs);
xc_std_lifscaled = std(xcs_lifscaled);
xcs_properr_m = mean(xcs_properr);
xcs_properr_std = std(xcs_properr);
xcs_perr_scaled_m = mean(xcs_properr_scaled);
xcs_perr_scaled_std = std(xcs_properr_scaled);
ac_ch1_mn = mean(ac_ch1);
ac_ch2_mn = mean(ac_ch2);
ac_ch1_std = std(ac_ch1);
ac_ch2_std = std(ac_ch2);

save(savename, 'xcmean', 'xc_std', 'ac_ch1_mn', 'ac_ch2_mn', 'ac_ch1_std', ...
    'ac_ch2_std', 'ave1', 'ave2', 'xcsmean_lifscaled', 'xcs_properr_m', ...
    'xcs_properr_std', 'xcs_perr_scaled_m', 'xcs_perr_scaled_std');