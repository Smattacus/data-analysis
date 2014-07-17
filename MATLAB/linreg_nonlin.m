function [ein, eout, mw_xf] = linreg_nonlin(npoints, nexps)
train = -1 + 2 * rand(2, 1000, nexps);
trainvals = sign(train(1,:,:).^2 + train(2,:,:).^2 -0.6);
%Implement noise:
trainvals = trainvals + (3 * rand(1, 1000, nexps) > .1);
trainvals = trainvals - (trainvals == 2) - 5 * (trainvals == 4);
%%%%
w_naive = zeros(3, nexps);
w_xf = zeros(6, nexps);
eins = zeros(1, nexps);
Y = zeros(1, 1000);
for i=1:nexps
    X = ones(3, 1000);
    X(2,:) = train(1, :, i);
    X(3,:) = train(2, :, i);
    Y = trainvals(1,:, i);
    w_naive(:, i) = ((X * X.')\X) * Y.';
    eins(i) = 1 / 1000 * sum((X.' * w_naive(:, i) - Y.').^2);
    X_xf = ones(6, 1000);
    X_xf(2,:) = train(1, :, i);
    X_xf(3,:) = train(2, :, i);
    X_xf(4,:) = train(1,:,i) .* train(2,:,i);
    X_xf(5,:) = train(1,:,i).^2;
    X_xf(6,:) = train(2,:,i).^2;
    w_xf(:, i) = ((X_xf * X_xf.')\X_xf) * Y.';
end
ein = mean(eins)
eouts = zeros(1, nexps);
mw_xf = mean(w_xf.');
for i=1:nexps
    tpts = (-1 + 2 * rand(2, 1000));
    Yo = sign(tpts(1, :).^2 + tpts(2,:).^2 - 0.6);
    Yo = Yo + 3 * (rand(1, 1000) > .1);
    Yo = Yo - (Yo == 2) - 5 * (Yo == 4);
    Xo_xf = ones(6, 1000);
    Xo_xf(2,:) = tpts(1,:);
    Xo_xf(3,:) = tpts(2,:);
    Xo_xf(4,:) = tpts(1,:) .* tpts(2,:);
    Xo_xf(5,:) = tpts(1,:).^2;
    Xo_xf(6,:) = tpts(2,:).^2;
    eouts(i) = (1 / 1000) * sum((Xo_xf.' * w_xf(:, i) - Yo.').^2);
end
mw_xf
eout = mean(eouts);
eout