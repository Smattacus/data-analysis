%%
%Generate the spectra for the surveys.
list = dir('*0,0--0,000*.txt');
for i=1:size(list,1)
    data = getLMPArray(list(i).name);
    AP1 = data{2};
    AP2 = data{3};
    xc1 = xcorr(AP1, AP2, 'unbiased');
    N = size(AP1, 1);
    t = (-(N-1):(N-1))/100000;
    [f, gx] = spec(xc1 .* exp(-(t/tc).^2/2), 1e-5);
    figure(1);
    semilogy(f, abs(gx));
    title(sprintf('Full spectrum for file %s', list(i).name);
    xlabel('Frequency');
    ylabel('|P|');
end