function [TOP1, BOT1, TOP2, BOT2] = genAllTPBox(path, phase, sq_freq, acq_f)
%Function to generate ALL the downsampled TOP and BOT arrays of LIF + Noise
%and noise data. 
%
%[TOP1, BOT1, TOP2, BOT2] = genAllTPBox(path, phase, freq, acq_f)
%
%This function assumes that the same phase is to be used across all files
%in the directory pointed to by path. It generates a straight square wave
%with phase and frequency given by phase and freq, respectively. This
%square wave is used to demodulate the data.
%
%This function finds all available files by searching for .h5 files. It
%will use all of them in the directory, so delete or move unwanted ones.
%
%INPUTS:
%   path - Path to a directory of h5 files.
%   phase - Phase to use in square wave function. 2x1 array for [pmt1,
%   pmt2].
%   freq - Frequency of the square wave function. Generally, this is the
%            laser chop frequency
%   acq_f   - Acquisition frequency.
%
%OUTPUTS:
%   TOP1 - NxM array of LIF +plasma data. N = number of files, M = # of points, PMT 1.
%   BOT1 - NxM array of plasma data, PMT 1.
%   TOP2 - NxM array of LIF + plasma data, PMT 2.
%   BOT2 - NxM array of plasma data, PMT 2.
cd(path);
list = dir('*.h5');

nfiles = size(list,1);
temp = h5read(list(1).name, '/PMT_DATA_8BIT');
N = size(temp,2);
total_t = N / acq_f;

base_phase = linspace(0, 2 * pi * total_t * sq_freq, N);
square1 = square(phase(1) + base_phase);
square2 = square(phase(2) + base_phase);

TOP1 = zeros(nfiles, N * sq_freq / acq_f);
BOT1 = zeros(nfiles, N * sq_freq / acq_f);
TOP2 = zeros(nfiles, N * sq_freq / acq_f);
BOT2 = zeros(nfiles, N * sq_freq / acq_f);

for i=1:size(list,1)
    fn = list(i).name;
    data = h5read(fn, '/PMT_DATA_8BIT');
    s1 = sum(data(1:16, :)); s1 = s1 - mean(s1);
    s2 = sum(data(17:32,:)); s2 = s2 - mean(s2);
    [T1, B1] = getTopBot(s1, square1, 1/acq_f, sq_freq/2);
    [T2, B2] = getTopBot(s2, square2, 1/acq_f, sq_freq/2);
    TOP1(i,:) = T1;
    BOT1(i,:) = B1;
    TOP2(i,:) = T2;
    BOT2(i,:) = B2;
end
    