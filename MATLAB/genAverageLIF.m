function genAverageLIF(path, savename, sq_freq, acq_freq)
%This function goes to a directory, pulls all the files, and generates the 
%average LIF photon count for each one. Then it writes the resulting array
%of values to the hard drive under 'savename'.
%INPUTS:
%   path - Path to the directory of the desired files.
%   savename - Name to save the files in.
%   sq_freq     - Frequency of square wave
%   acq_freq    - Frequency of acquisition frequency
cd(path)
[T1, B1, T2, B2] = genAllTPBoxFindPhase(path, sq_freq, acq_freq);
list = dir(sprintf('%s/*.h5', path));
D1 = T1 - B1;
D2 = T2 - B2;
ave1 = mean(D1.');
ave2 = mean(D2.');
save(savename, 'list', 'ave1', 'ave2');