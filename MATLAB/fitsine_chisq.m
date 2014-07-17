function [est, model] = fitsine_chisq(xdata, ydata, ydata_errors, start)
%Uses the fminsearch to fit a sine model to the data given by xdata and
%ydata. This implementation minimizes the chi squared function, rather than
%just the raw residuals.
%
%model =  A * sin ( 2 * pi * k * x + phi) + C
%
%function [est, model] = fitsine_chisq(xdata, ydata, ydata_errors, start)
%
%Inputs:
%xdata      - Vector of x coordinates for the data to be fit.
%ydata      - vector of y coordinates for the data to be fit.
%ydata_errors - Vector of errors corresponding to the y - values.
%start      - A 3 element vector containing a "good guess" for starting
%               parameters. [A, k, phi, C].
%Outputs:
%est        - Parameters which correspond to a local maximum found near the
%               starting point given by start.
%model      - The fitting function used. This is to make it easy to
%               generate a curve and residuals.
%
model = @fitfunc;
est = fminsearch(model, start);
    function [chi_sq, FittedCurve] = fitfunc(params)
        A = params(1);
        k = params(2);
        phi = params(3);
        C = params(4);
        FittedCurve = A * sin(2 * pi *  k * xdata + phi) + C;
        chi_sq = sum((1 ./ ydata_errors.^2).*(ydata - FittedCurve).^2);
    end
end