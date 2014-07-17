function gerrs = spec_errs(serrs)
%%
%Function to calculate the propagation of errors according to F. Skiff's
%spec function.
%
%gerrs = spec_errs(serrs, dt)
%
%Where serrs is the error array for the signal that was transformed. Note
%the signal array itself is not actually needed.
%
%
N = size(serrs,1);
if N == 1
    N = size(serrs,2);
end
%Use the downloaded function linspaceNDim to generate indices along
%each dimension to save on computation time
temp = linspace(1,N,N);
gerrs = zeros(1, N);
e = exp(-4i * pi / N * (temp - 1));
for i=1:N
    temp = sum(serrs.^2 .* e.^(i-1));
    gerrs(i) = 1 / sqrt(N) * sqrt(temp);
    if i < 20
        abs(gerrs(i))
    end
end