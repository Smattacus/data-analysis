function [TOP1, BOT1, TOP2, BOT2] = genAllTB(path, info_file, parallel, nworkers)
%Function to generate all of the downsampled TOP and BOTTOM matrices 
%of the LIF data.
%
%function [TOP1, BOT1, TOP2, BOT2] = genAllTB(path, parallel, nworkers, info_file)
%
%INPUTS:
% path          - Path to the desired folder.
% info_file     - Name of the file containing the center, delta info.
%                   Defaults to the first file with the suffix
%                   'chopinfo.txt' in the directory.
% parallel      - Paralelization enabled?
% nworkers      - Number of workers to open.
%
%OUTPUTS:
%TOPs           - NxM array of LIF ON data. N = number of files,
%                    M = downsampled length of array.
%BOTs           - NxM in the same way as TOPs. LIF OFF data.
%
cd(path);

chop_freq = 100000;
nyquist = 50000;
sample_dt = 1e-6;

if nargin < 3
    parallel = false;
    nworkers = 1;
end

if (nargin < 2)
    info_file = dir('*chopinfo.txt');
    info_file = info_file(1).name;
end

if parallel == true
    matlabpool('open', nworkers);
end

[fns centers deltas] = getChopInfo(info_file);
N = size(fns,1);
TOP1 = [];
BOT1 = [];
TOP2 = [];
BOT2 = [];
%For each file, demodulate the data and load the on and off data into the
%matrices.
for i=1:N
    file = fns{i};
    A = h5read(file, '/PMT_DATA_8BIT');
    s1 = sum(A(1:16,:));
    s1 = s1 - mean(s1);
    s2 = sum(A(17:32,:));
    s2 = s2 - mean(s2);
    sq1 = square(pi / 2 + genPhase(s1, centers(i), deltas(i), sample_dt));
    sq2 = square(pi / 2 + genPhase(s2, centers(i), deltas(i), sample_dt));
    [t1 b1] = getTopBot(s1, sq1, sample_dt, nyquist);
    [t2 b2] = getTopBot(s2, sq2, sample_dt, nyquist);
    TOP1 = [TOP1; t1];
    BOT1 = [BOT1; b1];
    TOP2 = [TOP2; t2];
    BOT2 = [BOT2; b2];
end
    
if parallel == true
    matlabpool('close');
end