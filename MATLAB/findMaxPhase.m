function [phase1, phase2] = findMaxPhase(filename, Ta, Fs, Fc, plots)
%Function which calculates and plots the maximum phase for the square wave.
%This assumes that you are using the synchronization box, and generates a
%square wave programmatically rather than using inverse transforms.
%
%It will display plots of the phase for the user to verify that the phase
%is correct and the same for both lasers.
%
%[phase1, phase2] = findMaxPhase(filename, Ta, Fs, Fc, plots)
%
%INPUTS:
%   filename
%   Ta          - Acq time in seconds
%   Fs          - Sampling frequency
%   Fc          - Chop frequency
%   plots       - Boolean to display plots
%
d = h5read(filename, '/PMT_DATA_8BIT');
s1 = sum(d(7:9,:));
s2 = sum(d(23:25,:));
phases = linspace(0, 2 * pi * Fc * Ta, Ta * Fs);
p = linspace(0, 2 * pi , 100);
PMT1 = zeros(1, 100);
PMT2 = zeros(1, 100);
for i=1:100
    PMT1(i) = sum((square(phases + p(i)) + 1) .* s1) / 2;
    PMT2(i) = sum((square(phases + p(i)) + 1) .* s2) / 2;
end
if plots == true
    figure(1); plot(p, PMT1); xlabel('Phase (rad)'); ylabel('a.u.'); title('PMT 1 Phase plot');
    figure(2); plot(p, PMT2); xlabel('Phase (rad)'); ylabel('a.u.'); title('PMT 2 Phase plot');
end
im1 = find(PMT1 == max(PMT1));
phase1 = p(im1);
im2 = find(PMT2 == max(PMT2));
phase2 = p(im2);