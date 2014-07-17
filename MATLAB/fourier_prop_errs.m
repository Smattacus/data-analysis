function [fbars] = fourier_prop_errs(fxn, means, errors, win, dt)
%This is a function which performs numerical error propagation in the
%manner outlined in Bevington & Robinson, 2nd ed, p 49-50. 
%
%This is designed for a Fourier transform, in which every data point is
%considered an independent source of error. A loop is created
%to iterate through the entire array array, adding error to one and only
%one point for each iteration of the loop.
%
%The resulting matrix of transforms is then summed in quadrature along 
%the relevant dimension and rooted to give an array of errors.
%
%[fbars] = fourier_prop_errs(fxn, means, errors, win, dt)
%
%INPUTS
%   fxn     - Function handle to some transform. FFT or spec recommended.
%   means   - Mean array.
%   errors  - The error bar array. 
%   win     - Windowing function, if desired. Set to 0 to ignore.
%   dt      - Time between samples.
if size(errors,1) == 1
    N = size(errors,2);
else
    N = size(errors,1);
end

props = zeros(N, 1).';

for i=1:N
    temp = means;
    temp(i) = temp(i) + errors(i);
    if win == 0
        [f, temp_arr] = fxn(temp, dt);
    else
        [f, temp_arr] = fxn(temp .* win, dt);
    end
    props = props + temp_arr.^2;
end
fbars = sqrt(props); %Sum over each column.