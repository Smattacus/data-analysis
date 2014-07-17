%%
%Quick routine to generate the power spectrum error bars according to the
%following process:
%   1) Load the means from the saved data from the Helium run.
%   2) Find the average of the standard deviation on the xcmean array.
%   3) Scale this average mean according to the ratio of window to box size
%   4) Scale this average according to the sqrt(n) of files in xcmean - 50
%   5) Plot this vs. the power spectrum of the data. 
%
%TODO:
%   1) Rerun some Helium runs while performing an error propagation
%   calculation to check the error arrays.
%   2) Do the above to generate the power spectrum error bars using the
%   theoretical errors from propagation and average LIF count.
%   3) Compare how these errors match up to the noise in the spectrum.
%   4) Repeat 1-3 for the normalized LIF data.
%
%%
list = dir('Dye_*.mat');
Nf = size(list,1);
N = 800000;
tc = 0.01;
xct = (-(N-1):(N-1))/1e5;
win = exp(-(xct/tc).^2/2);
winred = sqrt(sum(win.^2) / (2 * N - 1));
temp = ones(1, 2 * N - 1);
xcs_specs = zeros(21, 2 * N - 1);
xcs_rederrs = zeros(21,1);
%%
for i=1:Nf
    %%
    load(list(i).name);
    [f, gw] = spec(xcsmean_lifscaled .* win,1e-5);
    xcs_specs(i,:) = gw;
    e = mean(xc_std_lifscaled);
    e_red = winred * e / sqrt(size(ave1,2));
    xcs_rederrs(i) = e_red;
    figure(1);semilogy(f, abs(gw)); hold;
    semilogy(f, e_red * temp, 'r'); xlabel('Frequency'); ylabel('log(|P|)');
    legend('Power Spectrum', 'White noise level');
    title(sprintf('LIF Scaled Spectrum of %s', list(i).name));
    hold;
end
save('Specs_RedErrs.mat', 'list', 'xcs_specs', 'xcs_rederrs', 'f');