function [PMTdata, H] = getPMTData(filename, bool_8b)
%Function to readout all 32 channels of PMT data.
% [PMTdata, H] = getPMTData(filename, is8bit)
%Returns a 32xN data matrix of the PMT channels.
%Inputs:
%   filename - num0 file of data set.
%   is8bit  - Boolean whether the data is in 8bit format. True = 8bit
%   values extracted, while false = 32bit values extracted.
%Outputs:
%   PMTdata - Raw 32xN data
%   H       - Header from the data file.
if bool_8b == false
    [A1, H] = hread_PMTdata(filename);
elseif bool_8b == true
        [A1, H] = hread8bData(filename);
end
length = size(A1, 1);
acqs = length / 32;
PMTdata = reshape(A1.', 32, acqs);