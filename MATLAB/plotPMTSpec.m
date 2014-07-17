function [S, F, T, P] = plotPMTSpec(pmt_filename, hdf5)
%Function which plots the spectrograms of a given PMT file.
%
%[S, F, T, P] = plotPMTSpec(pmt_filename, hdf5)
%
%INPUTS:
%   pmt_filename - Filename of the relevant PMT file.
%   hdf5         - Boolean indicating whether the file is binary (false) or
%                   h5 format (true)
%OUTPUTS:
%   [S, F, T, P] - Standard outputs from spectrogram with 50 % overlap, 1/8
%   second spacing.
% 
%
    if hdf5 == false
        [A, H] = getPMTData(pmt_filename, true);
    elseif hdf5 == true
        A = h5read(pmt_filename, '/PMT_DATA_8BIT');
    end
    pmt1 = A(1:16, :);
    pmt2 = A(17:32,:);
    s1 = sum(pmt1); s1 = s1 - mean(s1);
    s2 = sum(pmt2); s2 = s2 - mean(s2);
    sums = [s1; s2];
    %Obtain a spectrogram with 50 % overlap and 1/8 second spacing:
    for j=1:2
        [S, F, T, P] = spectrogram(sums(j, :), 1048576/8, 524288/8, 1048576/8, 1000000);
        index = find(abs(F - 100000) > 400);
        Fsmall = F; Fsmall(index) = [];
        Psmall = ones(size(Fsmall,1), size(P,2));
        for i=1:size(P,2)
            temp = P(:,i); temp(index) = [];
            Psmall(:,i) = temp;
        end
        figure(j); surf(T, Fsmall, 10 * log(Psmall), 'edgecolor', 'none');
        xlabel('Time'); ylabel('Frequency'); 
        fn = strrep(pmt_filename, '_', '-');
        title(sprintf('Abbreviated Spectrogram of Power of PMT#%d\n File = %s', j, fn));
        view(0, 90)
    end
end