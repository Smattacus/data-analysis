function saveXCMeanBoxFindPhaseFastAcq(path, savename, sq_freq, acq_freq, PMTNum)
%This function generates and saves the averaged cross laser correlation
%function for a given folder.
%
%This is for fast acquisition data taken on 8 channels only. The PMT argument
%configures whether the 8 channels are just one PMT or two PMTs for purposes of
%phase calculation. In either case, channels 1-4 are cross correlated with
%channels 5-8.
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
% PMTNum        - Which PMT to use. Must take value of 1 or 2.
%
%
cd(path);
%[T1, B1, T2, B2] = genAllTB(path, info_file, parallel, nworkers);
%For now, use defaults:
[T1, B1, T2, B2] = genAllTPBoxFindPhaseFastAcq(path, sq_freq, acq_freq, 1);
list = dir(sprintf('%s/*.h5', path));
if PMTNum == 1
    D = T1 - B1;
elseif PMTNum == 2
    D = T2 - B2;
end
ave = mean(D.');
remove_ind = [];
%Remove data files that potentially have low LIF counts.
for i=1:size(D,1)
    if ave(i) < 0.1 
        display(sprintf('File %s omitted due to low LIF count.', list(i).name));
        display(sprintf('PMT1 = %d\nPMT2 = %d', ave(i)));
        remove_ind = [remove_ind, i];
    end
end
D(remove_ind,:) = [];
ave(remove_ind) = [];
N = size(D,1);
ac = zeros(size(D, 1), size(D,2) * 2 - 1);
ac_lifscaled = ac;
for i=1:size(D, 1);
   D(i,:) = D(i,:) - mean(D(i,:));
   ac(i,:) = xcorr(D(i,:), 'unbiased');
   ac_lifscaled(i,:) = xcorr(D(i,:) / ave(i), 'unbiased');
end
acmean = (mean(ac));
acmean_lifscaled = mean(ac_lifscaled);
ac_std = std(ac);
ac_std_lifscaled = std(ac_lifscaled);

save(savename, 'acmean', 'acmean_lifscaled', 'ac_std', 'ac_std_lifscaled');
