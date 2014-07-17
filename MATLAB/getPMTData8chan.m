function [PMTdata, H] = getPMTData8chan(filename)
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
[A1, H] = hread8bData(filename);
length = size(A1, 1);
acqs = length / 8;
PMTdata = reshape(A1.', 8, acqs);
