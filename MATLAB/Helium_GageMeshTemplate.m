%%
%This just does the averaged cross correlation across 40 files of channel B
%of the gage scope. It plots this graph to file.
%It will also make a mesh plot of the individual gage scope spectra from 0
%- 7500 Hz.
%%
%Set up some file stuff.
addpath('/Users/smattingly/Programs/MATLAB');
dye = 'dyenamehere';
datapath = '/Users/smattingly/Data/11-5-2013/V2PT';
gagemesh = sprintf('%s/GageMeshes', datapath);
path = sprintf('%s/%s', datapath, dye);
path
cd(path)
mkdir GageMeshes
%%
list = dir(sprintf('*%s*.h5', strrep(dye, '.', 'p')));
temp = h5read(list(1).name, '/GAGE_CHAN_B');
gagecorr = zeros(size(list,1), 2 * size(temp,1) - 1);
gagespec = gagecorr;
N = size(temp,1);
t = (-(N-1):(N-1)).'/1e5;
Nf = size(list,1);
first = list(1).name;
present = false;
for i=1:Nf
    %Make logic to use the available files, not just construct the names
    %blindly
    fn = strrep(first, 'num0', sprintf('num%d', i-1));
    present = false;
%
%This logic is to prevent deleted files from screwing up the loop.
    for j=1:Nf
        if strcmp(fn, list(j).name);
            present = true;
            break;
        end
    end
    if present == false
        continue;
    end
%End of logic.
    a = h5read(fn, '/GAGE_CHAN_B');
    ac = xcorr(double(a - mean(a)), 'unbiased');
    gagecorr(i, :) = ac;
    [f, gagespec(i,:)] = spec(ac .* exp(-(t/.01).^2/2), 1e-5);
end
%%
gcmn = mean(gagecorr, 1).';
it = find(abs(t) < .05);
[f, g] = spec(gcmn .* exp(-(t/.01).^2/2), 1e-5);
ig = find(abs(f) < 12500); 
figure(2); semilogy(f(ig), abs(g(ig))); xlabel('Frequency'); ylabel('log abs power');
title('Spectrum of averaged auto correlation for 10 W, SP780, 200utorr');

h1 = figure(1); set(h1, 'visible', 'off');
semilogy(f, abs(g)); title(sprintf('Averaged Gagescope Spectrum for\nPath = %s', path));
xlabel('f (Hz)');
ylabel('Abs power');
xcfn = sprintf('GageMeshes/AVERAGED_GAGE_CHANB.png');
print(h1, '-dpng', xcfn);

h2 = figure(2); set(h2, 'visible', 'off');
semilogy(f(ig), abs(g(ig))); title(sprintf('Averaged Gagescope Spectrum < 7500 Hz for path =%s', path));
xlabel('f (Hz)');
ylabel('Abs power');
xcfn = 'GageMeshes/AVERAGED_GAGE_CHANB_ZOOM.png';
print(h2, '-dpng', xcfn);



%%
%Make a mesh of the spectra.
h3 = figure(3); set(h3, 'visible', 'off');
ig = find(abs(f) < 7500); 
mesh(f(ig), 1:Nf, abs(gagespec(:, ig))); view(0,90); xlabel('Freq (Hz)'); ylabel('Acq #');
title(sprintf('Spectrogram of gage acqs for path = %s', path));
xcfn = 'GageMeshes/GAGE_SPECTROGRAM.png';
print(h3, '-dpng', xcfn);
%%
exit;