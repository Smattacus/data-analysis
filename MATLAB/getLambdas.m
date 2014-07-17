function [dye, diode] = getLambdas(fn)
%Function for deriving the dye and diode values from the filename.
%
% [dye_wavelength, diode_wavelength] = getLambdas(fn)
%
% This code assumes that the filenames are of the format
%
%Mon_Aug_05_21_16_15_num0_ZeemanScan_11W_150utorr_Dye_611-66780_Diode_668-61020.h5 
%
% Where the wavelength values are filled in as "Dye_xxx-yyyyy" and
% "Diode_xxx-yyyyy."
%
i_dye = strfind(fn, 'Dye_') + 4;
dye = str2num(sprintf('%s.%s', fn(i_dye:i_dye+2), fn(i_dye+4:i_dye+8)));
i_diode = strfind(fn, 'Diode_') + 6;
diode = str2num(sprintf('%s.%s', fn(i_diode:i_diode+2), fn(i_diode+4:i_diode+8)));