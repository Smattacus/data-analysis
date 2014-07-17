function [ein, eout] = linreg_test(npoints, nexps)
%Script for implementing the linear regression for CS156 HW#2.
%Define a target function (Evaluate it a thousand times!):
%x1, y1, x2, y2
rand_targets = -1 + 2 * rand(4, nexps);
slopes = (rand_targets(4, :) - rand_targets(2,:)) ./ ...
        (rand_targets(3,:) - rand_targets(1,:));
ints = rand_targets(2, :)  - slopes .* rand_targets(1,:);
ein = zeros(nexps,1);
eout = zeros(nexps,1);
%Generate random data set.
%2 coords, npoints data points, nexps sets of data.
for i=1:nexps
    %(x1, y1), npoints columns
    r = -1 + 2 * rand(2, npoints);
%Plug into their respective target functions
    dataset = ((ints(i) + slopes(i) .* r(1, :)) - r(2, :) ) < 0;
    %Assign zeros to -1
    dataset = dataset + (dataset == 0) * -1;
    %Using the linreg algorithm from class. Our definitions are already
    %transposed:
    w = ((r * r.')\ r) * dataset.';
    %Fitted set
    fitset = ((r.' * w).' < 0);
    fitset = fitset + (fitset == 0) *-1;
    ein(i) = (1 / npoints) * sum(fitset == dataset);
%    fitset(i,1) = 1 / npoints * (r.' * w - dataset.').' * (r.' * w - dataset.');
    %Problem 6:
    %Generate 1000 random points:
    rout = -1 + 2 * rand(2, npoints);
    tout = ((ints(i) + slopes(i) .* rout(1,:)) - r(2,:)) < 0;
    tout = tout + (tout == 0) * -1;
    fitout = ((rout.' * w) < 0).';
    fitout = fitout + (fitout == 0 ) * -1;
    eout(i) = 1 / npoints * sum(fitout == tout);
end
ein = mean(ein);
eout = mean(eout);

