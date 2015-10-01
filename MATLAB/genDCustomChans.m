function LIFS = genDCustomChans(fname, Ta, Fs, Fc, varargin)
%Function which generates the average LIF photon rate per demodulated sample
%time for the sets of PMT channels as given by args.
%
%LIFS = genDCustomChans(fname, Ta, Fs, Fc, args)
%
%INPUTS:
%   fname   :   Filename of the desired data file.
%   Ta      :   Total time (in seconds) of acquired data.
%   Fs      :   Sampling frequency.
%   Fc      :   Chop frequency.
%   varargin:   Variable arguments to specify PMT channels. See below.
%
%OUTPUTS:
%   LIFS    :   Array of average LIF photon rate per demodulated sample (per
%   chop period)
%

%Number of variable inputs
Nv = - nargin - 4;

%Load the data and find the phase.
data = h5read(fname, '/PMT_DATA_8BIT');
N = size(data,2);
total_t = N/Fs;

s = zeros(size(data,2), 1)

base_phase = genBasePhase(total_t);

[p1, p2] = findMaxPhase(fname, total_t, Fs, Fc, true);
%This assumes that size(p1) = size(p2) = 1. If it's not, the file is probably
%borked.
sq1 = square(p1 + base_phase);
sq2 = square(p2 + base_phase);
display(sprint('Phase 1 = %f', p1));
display(sprint('Phase 2 = %f', p2));

diffs = zeros(Nv, 1)

%Iterate through each entry of varargin.
for i=1:Nv
    chans = varargin(i);
    s = zeros(size(data,2), 1);
    if ischar(chans)
        %Assumed to be of form 'N1-N2'
        %e.g. 8-15
        Cb = strsplit(chans, '-');
        N1 = Cb(1);
        N2 = Cb(2);
        s = sum(data(N1:N2, :));
        %WARNING: Using phase generated from PMT #1!!!.
        if N1 <= 16 && N2 <= 16
            [T, B] = getTopBot(s, sq1, 1/Fs, Fc/2);
        else if N1 > 16 && N2 > 16
            [T, B] = getTopBot(s, sq2, 1/Fs, Fc/2);
        else
            %If N1 and N2 ranges over two difference PMTs, then just
            %use the square wave for PMT1 (hopefully they're the same).
            [T, B] = getTopBot(s, sq1, 1/Fs, Fc/2);
            display('Beware! Using phase from PMT1 for channel sum across PMT1 & PMT 2')
            D = T - B;
            diffs(i) = mean(D);
        end
    else
        cs = size(chans);
        cs = cs(1);
        for j=1:cs
            s = s + data(chans(j), :);
        end
        if chans(1) <= 16
            [T, B] = getTopBot(s, sq1, 1/Fs, Fc/2);
        else
            [T, B] = getTopBot(s, sq2, 1/Fs, Fc/2);
        D = T - B;
        diffs(i) = mean(D);
    end
end
return diffs


