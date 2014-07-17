%Script for generating all of the PMT Spectrograms in order to determine
%good windows for the chop frequency (or if a file should be rejected since
%the LIF signal is too weak.
%%
list = dir('*.DAT');
for i=1:size(list,1);
%%
    fn = list(i).name;
    plotPMTSpec(fn);
end

%Analyze the data, file by file:
%%
%Generate the square wave demodulations.
S1 = zeros(12, size(list,1));
S2 = zeros(12, size(list,1));
for i=1:size(list,1);
    fn = list(i).name;
    [A, H] = getPMTData(fn, true);
    s1 = sum(A(1:16,:));
    s2 = sum(A(17:32,:));
    ph1 = genPhase(s1, 1025, 30, 1e-6);
    ph2 = genPhase(s1, 1025, 30, 1e-6);
    %Check which phase gives a maximum for the demodulated signal
    for j=1:12
        S1(j, i) = sum(square(ph1 + j * pi * 2 / 12) .* s1);
        S2(j, i) = sum(square(ph1 + j * pi * 2 / 12) .* s2);
    end
%%
figure(1); plot(S1(:,i)); title('PMT 1 square wave corr');
figure(2); plot(S2(:,i)); title('PMT 2 square wave corr');
end
%The phase required is pi / 2 for all the data. Let's do it and generate
%xcorrs:
%%
clear S1 S2
tc = 0.01;
%Adjust for file size
tp = zeros(size(list,1), 3276799);
for i=1:size(list,1)
    ftitle = strrep(list(i).name, '_', '-');
    [A,H] = getPMTData(list(i).name, true);
    s1 = sum(A(1:16,:));
    s2 = sum(A(17:32,:));
    sq1 = square(pi / 2 + genPhase(s1, 100010, 50, 1e-6));
    sq2 = square(pi / 2 + genPhase(s2, 100010, 50, 1e-6));
    [TOP1, BOT1] = getTopBot(s1, sq1, 1e-6, 50000);
    [TOP2, BOT2] = getTopBot(s2, sq2, 1e-6, 50000);
    DIFF1 = TOP1 - BOT1; 
    DIFF2 = TOP2 - BOT2;
    N = size(TOP1, 2);
    t = (-(N-1):(N-1))/100000;
    %Windowing function
    win = exp(-(t/tc).^2/2);
%%
    %Autocorrelations of demodulated data to check
    x1 = xcorr(TOP1, 'unbiased');
    x2 = xcorr(BOT1, 'unbiased');
    x3 = xcorr(DIFF1, 'unbiased');
    [f, gT1] = spec(x1 .* win, 1e-5);
    [f, gB1] = spec(x2 .* win, 1e-5);
    [f, gD1] = spec(x3 .* win, 1e-5);
    ig = find(abs(f) < 5000);
    figure(1); semilogy(f(ig), abs(gT1(ig))); title('ACORR PMT 1 TOP');
    figure(2); semilogy(f(ig), abs(gB1(ig))); title('ACORR PMT 1 BOT');
    figure(3); semilogy(f(ig), abs(gD1(ig))); title('ACORR PMT 1 DIFF');
%%
    x1 = xcorr(TOP2, 'unbiased');
    x2 = xcorr(BOT2, 'unbiased');
    x3 = xcorr(DIFF2, 'unbiased');
    [f, gT2] = spec(x1 .* win, 1e-5);
    [f, gB2] = spec(x2 .* win, 1e-5);
    [f, gD2] = spec(x3 .* win, 1e-5);
    figure(1); semilogy(f(ig), abs(gT2(ig))); title('ACORR PMT 2 TOP');
    figure(2); semilogy(f(ig), abs(gB2(ig))); title('ACORR PMT 2 BOT');
    figure(3); semilogy(f(ig), abs(gD2(ig))); title('ACORR PMT 2 DIFF');
%%
    %Cross correlation
    xD = xcorr(DIFF1, DIFF2, 'unbiased');
    [f, gX] = spec(xD .* win, 1e-5);
    figure(1); semilogy(f(ig), abs(gX(ig)));
    title(sprintf('2 Point Correlation Spectrum\nFor file %s', ftitle));
    %Save the cross correlation to an array for averaging later.
    tp(i,:);
end    