%Script to generate mesh and contour plots to files for a list of files
%in a given text file.
%For now the input file is all_arrays.txt

f = fopen('all_arrays.txt');
tline = strtrim(fgetl(f));
freq = 55.5e3

while tline ~= -1
    [fr, Q, H, data] = lmpSpec2(tline, .015, freq);
    LP = data{:,1};
    AP1 = data{:,2};
    AP2 = data{:,3};
    c1 = xcorr(LP, 'unbiased');
    c2 = xcorr(LP, AP1, 'unbiased');
    s = size(LP, 1);
    t = (-(s(1) - 1):s(1) - 1) / freq; %Create the time axis.
    [f1, g1] = spec(c1 .* exp(-(t' / .015).^2/2), 1/freq);
    [f2, g2] = spec(c2 .* exp(-(t' / .015).^2/2), 1/freq);
    %Mesh plot
    display('Creating mesh plot')
    figure(1); 
    mesh(-4:3, fr(56002:end), Q(56002:end,:));
    display('Mesh plot created')
    xlabel('Mode number');
    ylabel('Frequency');
    zlabel('Log normal power');
    angle = getAngle(tline);
    suff= strrep(getSuffix(tline), '.txt', '.jpg');
    suff_es = strrep(suff, '_', '\_');
    tn = sprintf('Spectral Mesh plot for theta = %f deg, suffix = %s', angle, suff_es);    
    fn = sprintf('MeshTheta=%fSuffix=%s', angle, suff);
    title(tn);
    print(fn, '-djpeg');
    %Contuor plot
    contour(-4:3, fr(56002:end), Q(56002:end,:));
    colorbar
    xlabel('Mode number');
    ylabel('Frequency');
    angle = getAngle(tline);
    tn = sprintf('Spectral contuor plot for theta = %f deg, suffix = %s', angle, suff_es);    
    fn = sprintf('ContuorTheta=%fSuffix=%s', angle, suff);
    title(tn);
    print(fn, '-djpeg');
    %Dual plot of autocorrelat and cross correlated langmuir probe (LP).
    semilogy(f1, abs(g1).^2, f1, abs(g2).^2);
    xlabel('Frequency (Hz)');
    ylabel('Power squared spectrum');
    tn = sprintf('Power Spectrum for theta = %f deg, suffix = %s', angle, suff_es);    
    title(tn);
    legend('Probe 1 Auto', 'Probe 1, 2 cross');
    fn = sprintf('PowerSpecAng=%fSuffix=%s', angle, suff);
    print(fn, '-djpeg');
    %
    %Generate the next iteration's filename.
    %
    nline = fgetl(f);
    if nline == -1
        break
    end
    tline = strtrim(fgetl(f));
end