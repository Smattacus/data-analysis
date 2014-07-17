function [X V S U] = genSVDFromArray(filename)
%Function which takes in an LMP array + standalone file and
%and generates the rectangular matrix which can be decomposed with SVD.
%
% [X] = genSVDData(filename)
%
%INPUTS:
% FILENAME  - Text file of a LMP array data set.
% 
%OUTPUTS:
% X         - Square matrix of array data, with rows in time and columns in
% space (i.e. Nsamples x 8)
%
data = getLMPArray(filename);
%Create a matrix of Nsamp x 8 columns
X = [data{2}, data{3}, data{4}, data{5}, data{6}, data{7}, data{8}, data{9}];
%I am switching the V and U to match the notation in C Nardone's paper
%from 1992, Plasma Phys Control Fusion 34 1447.
[V S U] = svd(X, 'econ');