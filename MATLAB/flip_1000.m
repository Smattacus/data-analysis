function [v1, vrand, vmin] = flip_1000()
%Script to implement Hoeffding inequality code
flips = (rand([10, 1000]) > 0.5);
size(flips);
flipsums = sum(flips);
[minval, cmin] = min(flipsums);
%Assign the v's:
v1 = flipsums(1) / 10;
vrand = flipsums(randi(1000, 1, 1)) / 10;
vmin = minval / 10;