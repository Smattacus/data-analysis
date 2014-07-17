function datamatrix = createDM(filename, submean, tc)
% Quick function to create a data matrix for the desired .MAT file. Windows
%according to tc. The output datamatrix is of the form [f, abs(g)] for the 
%first 10000 Hz.
%
% datamatrix = createDM(filename, submean, tc)
%
% INPUTS
% filename      - Desired .MAT file of the xcmean
% submean       - Boolean. TRUE - subtracts the mean. FALSE - no change.
% tc            - Correlation time to window by.
% OUTPUTS
% datamatrix    - Variable to plug into IPF.
S = load(filename);
xcmean = S.xcmean;
if submean == true
    xcmean = xcmean - mean(xcmean);
end
N = (size(xcmean,2) + 1) / 2;
t = (-(N-1):(N-1))/1e5;
[f, g] = spec(xcmean .* exp(-(t/tc).^2/2), 1e-5);
ig = find(abs(f) < 10000);
datamatrix = [f(ig); abs(g(ig))];
end