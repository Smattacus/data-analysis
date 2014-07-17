function [n, Vp, Te] = read_params(filename)
%Function to read a LMP parameters text file.
%Inputs:
%   Filename    - Filename of the LMP param file. usually *_param.txt.
%Outputs:
%   n           - Number density of the plasma in cm^-3
%   Vp          - Plasma potential
%   Te         -  Electron temperature