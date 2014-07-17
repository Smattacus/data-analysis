function xc_errs = xcorr_err(sig1, sig2, sig1err, sig2err)
%Function to generate the analytical values from a error propagation
%formula applied to MATLAB's cross correlation function.
%
%This function returns an array of errors corresponding to the errors at
%each cross correlation point after a cross correlation.
%
%xc_errs = xcorr_err(sig1, sig2, sig1err, sig2err)
%
%This function assumes that the error on each point in the signals is
%constant across the signal array. More in depth analysis will require that
%this is not the case.
%

N = size(sig1,2);
if N ==1
    N = size(sig1,1);
end

xc_errs = zeros(1, 2 * N - 1);
%This is broken up into three parts for now.
%Center of array (n' = N-1)
var = 1 / N^2 * (sig1err^2 * sum(sig2.^2) + sig2err^2 * sum(sig1.^2));
%Sides of correlation
xc_errs(N) = sqrt(var);
for i=1:N-1
    temp1 = 1/(N-i)^2 * (sig1err^2 * sum(sig2(1:N-i).^2) +  ...
        sig2err^2 * sum(sig1(i:N).^2));
    xc_errs(N+i) = sqrt(temp1);
    temp2 = 1 / (N-i)^2 * (sig1err^2 * sum(sig2(i:N).^2) + ...
        sig2err^2 * sum(sig1(1:N-i).^2));
    xc_errs(N-i) = sqrt(temp2);
end
