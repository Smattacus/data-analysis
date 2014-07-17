function v = dopplerVelocity(lambda_0, lambda)
%Calculates the Doppler velocity according to the Doppler equation
%
%sigma = sigma_0  * ( 1 + v / c)
%
%Where sigma = 1 / wavelength. Enter lambda_0 and lambda in nanometers.
%
sigma = 1 ./ lambda;
sigma_0 = 1 / lambda_0;
v = 2.9979e8 * (sigma / sigma_0 - 1);