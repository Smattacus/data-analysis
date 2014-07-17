function printAllPhasesBox(path)
%This function generates the square wave corresponding to the chopped
%signal and plots the correlation over a single period.
%
% printAllPhases(path, parallel
%
% INPUTS: 
%   path     - Desired path to run in.
%   parallel - Boolean to determine whether to use parallel statements
%   nworkers - Number of workers to start up.
%
%This function assumes that the synchronization box has been used in order
%to make sure that the laser chop is at 100 khz. This function is best used
%to verify the phase of the laser in relation to a straight 100 khz square
%wave.
%
%

cd(path);
if exist('PhasePNG', 'dir') == 0
    mkdir PhasePNG
end
fns = dir('*.h5');

%[fns centers deltas] = getChopInfo(list(1).name);
numcorrs = 100;
S1 = zeros(size(fns,1), numcorrs+1);
S2 = zeros(size(fns,1), numcorrs+1);
phase = linspace(0, 2 * pi * 8 * 100e3, 1e6 * 8);
phaxis = linspace(0, 2 * pi, numcorrs+1);

for i=1:size(fns,1)
    fn = fns(i).name;
    A = h5read(fn, '/PMT_DATA_8BIT');
    s1 = sum(A(1:16, :));
    s1 = s1 - mean(s1);
    s2 = sum(A(17:32, :));
    s2 = s2 - mean(s2);
    %Check which phase gives a maximum for the demodulated signal
    for j=1:numcorrs+1
        [TOPS1, BOTS1] = getTopBot(s1, square(phase + 2 * pi * (j-1)/numcorrs), 1e-6, 50000);
        [TOPS2, BOTS2] = getTopBot(s2, square(phase + 2 * pi * (j-1)/numcorrs), 1e-6, 50000);        
        S1(i, j) = sum(TOPS1 - BOTS1);
        S2(i, j) = sum(TOPS2 - BOTS2);
    end
    h1 = figure(1); set(h1, 'visible', 'off');
    plot(phaxis, S1(i,:)); title('PMT #1 Square Wave Corr');
    xlabel(' Displacement (radians)');
    ylabel('Correlation');
    xcfn = sprintf('PhasePNG/%s_PHASE_PMT1.png', strrep(fn, '.h5', ''));
    legend(sprintf('Max at %f rad', phaxis(S1(i,:) == max(S1(i,:)))));
    print(h1, '-dpng', xcfn);
    h2 = figure(2); set(h2, 'visible', 'off');
    plot(phaxis, S2(i,:)); title('PMT #2 Square Wave Corr');
    xlabel('Displacement (2 pi / 12)');
    ylabel('Correlation');
    legend(sprintf('Max = %f', phaxis(S2(i,:) == max(S2(i,:)))));
    xcfn = sprintf('PhasePNG/%s_PHASE_PMT2.png', strrep(fn, '.h5', ''));
    print(h2, '-dpng', xcfn);
end

%{
if parallel == true
    matlabpool('open', nworkers);
    
    parfor i=1:size(list, 1)
        fn = fns[i];
        A = h5read(fn, '/PMT_DATA_8BIT');
        s1 = sum(A(1:16, :));
        s2 = sum(A(17:32, :));
        ph1 = genPhase(s1, centers(i), deltas(i), 1e-6);
        ph2 = genPhase(s2, centers(i), deltas(i), 1e-6);
        %Check which phase gives a maximum for the demodulated signal
        for j=1:12
            S1(i, j) = sum(square(ph1 + j * pi * 2 / 12) .* s1);
            S2(i, j) = sum(square(ph2 + j * pi * 2 / 12) .* s2);
        end
        h1 = figure; 
        h1 = plot(S1(i,:)); title('PMT #1 Square Wave Corr'); 
        xlabel(' Displacement (2pi / 12)');
        ylabel('Correlation');
        xcfn = sprintf('PhasePNG/%s_PHASE_PMT1.png', strrep(fn, '.h5', ''));
        print(h1, '-dpng', xcfn);
        h2 = figure;
        h2 = plot(S2(:,i); title('PMT #2 Square Wave Corr');
        xlabel('Displacement (2 pi / 12)');
        ylabel('Correlation');
        xcfn = sprintf('PhasePNG/%s_PHASE_PMT2.png', strrep(fn, '.h5', ''));
        print(h2, '-dpng', xcfn);
    end
    matlabpool('close');
    %}
