function [TOP1, BOT1, TOP2, BOT2] = genAllTPBoxFindPhase_Cut(path, Fc, Fa,
varargin)
%Function to generate ALL the downsampled TOP and BOT arrays of LIF + Noise
%and noise data. 
%
%This one is different in that it drops the first 300k elements. We're
%going to see if this helps to remove the transient behavior.
%
%[TOP1, BOT1, TOP2, BOT2] = genAllTPBoxFindPhase(path, Fc, Fa, start)
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
%   Fc - Frequency of the square wave function. Generally, this is the
%            laser chop frequency
%   Fa   - Acquisition frequency.
%   start   -   (OPTIONAL) Specifies the index to start at for the data arrays.
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
total_t = N / Fa;

if nVarargs < 1
    start = 300000
else nVarargs == 1
    start = varargin{1}
else 
    display('Incorrect number of arguments in genAllTPBoxFindPhase_Cut')
end

base_phase = genBasePhase(total_t);

TOP1 = zeros(nfiles, N * Fc / Fa);
BOT1 = zeros(nfiles, N * Fc / Fa);
TOP2 = zeros(nfiles, N * Fc / Fa);
BOT2 = zeros(nfiles, N * Fc / Fa);

for i=1:size(list,1)
    fn = list(i).name;
    [p1, p2] = findMaxPhase(fn, total_t, Fa, Fc, true);
    if size(p1,1) > 1 || size(p1,2) > 1
        %This should only happen if p1 or p2 is = 0.
        display(sprintf('Phase longer than one element for file = %s, Ch1.', fn));
        display('Taking first element of phase 1 array');
        p1 = p1(1);
    end
    if size(p2,1) > 1 || size(p2,2) > 1
        %This should only happen if p1 or p2 is = 0.
        display(sprintf('Phase longer than one element for file = %s, Ch1.', fn));
        display('Taking first element of phase 2 array');
        p2 = p2(1);
    end
    square1 = square(p1 + base_phase);
    square2 = square(p2 + base_phase);
    display(p1);
    display(p2);
    data = h5read(fn, '/PMT_DATA_8BIT');
    s1 = sum(data(1:16, start:)); 
%    s1 = s1 - mean(s1);
    s2 = sum(data(17:32,start:)); 
%    s2 = s2 - mean(s2);
    [T1, B1] = getTopBot(s1, square1, 1/Fa, Fc/2);
    [T2, B2] = getTopBot(s2, square2, 1/Fa, Fc/2);
    TOP1(i,:) = T1;
    BOT1(i,:) = B1;
    TOP2(i,:) = T2;
    BOT2(i,:) = B2;
end
