function [A, H] = hread_PMTdata(filename)
%Function to read the PMT data read with the VME acquisition system.
%
% [A, H] = hread_PMTdata(filename);
%
%Inputs:
%filename   - *.DAT output from the C function write_VME
%Outputs:
%A          - PMT data array.
%H          - Header elements describing acquisition time.
%

f1 = fopen(filename);
H = [fread(f1, 128, '*char', 'l').'];
for i=1:7
    H = [H; fread(f1, 128, '*char').'];
end
A = fread(f1, '*uint');
end