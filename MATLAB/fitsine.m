function [est, model] = fitsine(xdata, ydata, start)
%Uses the fminsearch to fit a sine model to the data given by xdata and
%ydata. This implementation does not use error bars. The sine function used
%is 
%
%model =  A * sin ( 2 * pi * k * x + phi) + C
%
%[est, model] = fitsine(xdata, ydata, start)
%
%Inputs:
%xdata      - Vector of x coordinates for the data to be fit.
%ydata      - vector of y coordinates for the data to be fit.
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
    function [sse, FittedCurve] = fitfunc(params)
        A = params(1);
        k = params(2);
        phi = params(3);
        C = params(4);
        FittedCurve = A * sin(2 * pi *  k * xdata + phi) + C;
        sse = sum((ydata - FittedCurve).^2);
    end
end