function saveXCMeanBoxFindPhase(path, savename, sq_freq, acq_freq, varargin)
%This function generates and saves the averaged cross laser correlation
%function for a given folder.
%
%This function also saves the mean LIF photon count per file in a separate
%array.
%
% saveXCMeanBoxFindPhase(path, savename, sq_freq, acq_freq)
%
% INPUTS:
% path          - Directory to save the *.mat file to.
% savename      - Will be used to create the .mat file, according to
%                   savename.mat
% sq_freq       - Frequency of laser chop
% acq_freq      - Frequency of acquisition
%start          - Location to eliminate data.
%
cd(path);
%[T1, B1, T2, B2] = genAllTB(path, info_file, parallel, nworkers);
%For now, use defaults:
nvargs = length(varargin)
[T1, B1, T2, B2] = genAllTPBoxFindPhase(path, sq_freq, acq_freq);

if nvargs == 1
    start = varargin{1};
    T1 = T1(:, start:end);
    T2 = T2(:, start:end);
    B1 = B1(:, start:end);
    B2 = B2(:, start:end);
end

list_data = dir(sprintf('%s/*.h5', path));
D1 = (T1 - B1) * 2;
D2 = (T2 - B2) * 2;
ave1 = mean(D1.');
ave2 = mean(D2.');
bg_ave1 = mean(B1.');
bg_ave2 = mean(B2.');
remove_ind = [];
%Remove data files that potentially have low LIF counts.
for i=1:size(D1,1)
    if ave1(i) < 0.01 || ave2(i) < 0.01
        display(sprintf('File %s omitted due to low LIF count.', list_data(i).name));
        display(sprintf('PMT1 = %d\nPMT2 = %d', ave1(i), ave2(i)));
        remove_ind = [remove_ind, i];
    end
end
D1(remove_ind,:) = [];
D2(remove_ind,:) = [];
ave1(remove_ind) = [];
ave2(remove_ind) = [];
list_data(remove_ind) = [];
N = size(D1,1);
xcs = zeros(size(D1, 1), size(D1,2) * 2 - 1);
xcs_lifscaled = xcs;
ac_ch1 = xcs;
ac_ch2 = xcs;
%Adding this in 3 Mar 2017 to estimate background correlations
%and subtract them. SWM.
xcs_background = xcs;
xc_bg1_lif2 = xcs;
xc_bg2_lif1 = xcs;
D1_scale = D1 * 0;
D2_scale = D2 * 0;
for i=1:size(D1, 1);
   D1(i,:) = D1(i,:) - mean(D1(i,:));
   D2(i,:) = D2(i,:) - mean(D2(i,:));
   xcs(i,:) = xcorr(D1(i,:), D2(i,:), 'unbiased');
   xcs_lifscaled(i,:) = xcorr(D1(i,:) / ave1(i), D2(i,:) / ave2(i), 'unbiased');
   %3 Mar 2017
   xcs_background(i,:) = xcorr((B1(i,:) - bg_ave1(i))/bg_ave1(i), (B2(i,:) - bg_ave2(i))/bg_ave2(i), 'unbiased');
   xc_bg1_lif2(i,:) = xcorr( (B1(i,:) - bg_ave1(i)) / bg_ave1(i), (D2(i,:) / ave2(i)), 'unbiased');
   xc_bg2_lif1(i,:) = xcorr( (D1(i,:) / ave1(i)), (B2(i,:) - bg_ave2(i)) / bg_ave2(i), 'unbiased');
   ac_ch1(i,:) = xcorr(D1(i,:), 'unbiased');
   ac_ch2(i,:) = xcorr(D2(i,:), 'unbiased');
   D1_scale(i,:) = D1(i,:) / ave1(i);
   D2_scale(i,:) = D2(i,:) / ave2(i);
end
xcmean = (mean(xcs));
xcsmean_lifscaled = mean(xcs_lifscaled);
xc_std = std(xcs);
xc_std_lifscaled = std(xcs_lifscaled);
ac_ch1_mn = mean(ac_ch1);
ac_ch2_mn = mean(ac_ch2);
ac_ch1_std = std(ac_ch1);
ac_ch2_std = std(ac_ch2);
xcm_bg12 = mean(xcs_background);
xcm_bg12_std = std(xcs_background);
xcm_bg1_lif2 = mean(xc_bg1_lif2);
xcm_bg1_lif2_std = std(xc_bg1_lif2);
xcm_bg2_lif1 = mean(xc_bg2_lif1);
xcm_bg2_lif1_std = std(xc_bg2_lif1);

save(savename, 'xcmean', 'xc_std', 'ac_ch1_mn', 'ac_ch2_mn', 'ac_ch1_std', ...
    'ac_ch2_std', 'ave1', 'ave2', 'xcsmean_lifscaled', 'xc_std_lifscaled', ...
    'xcm_bg12', 'xcm_bg12_std','xcm_bg1_lif2', 'xcm_bg1_lif2_std', 'xcm_bg2_lif1', 'xcm_bg2_lif1_std', ...
    'bg_ave1','bg_ave2','list_data');
%Secondary: Save the normalized D1 and D2 matrices.
savename_diffs = strrep(savename, '.mat', '_DiffArrays.mat');
save(savename_diffs, 'D1_scale', 'D2_scale', 'ave1', 'ave2');
