function printAllPhases(path, parallel, nworkers)
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
cd(path);
if exist('PhasePNG', 'dir') == 0
    mkdir PhasePNG
end
list = dir('*chopinfo.txt');

if nargin == 1
    parallel = false;
    nworkers = 1;
end

[fns centers deltas] = getChopInfo(list(1).name);
S1 = zeros(size(fns,1), 12);
S2 = zeros(size(fns,1), 12);

for i=1:size(fns,1)
    fn = fns{i};
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
    h1 = figure(1); set(h1, 'visible', 'off');
    plot(S1(i,:)); title('PMT #1 Square Wave Corr');
    xlabel(' Displacement (2pi / 12)');
    ylabel('Correlation');
    xcfn = sprintf('PhasePNG/%s_PHASE_PMT1.png', strrep(fn, '.h5', ''));
    print(h1, '-dpng', xcfn);
    h2 = figure(2); set(h2, 'visible', 'off');
    plot(S2(i,:)); title('PMT #2 Square Wave Corr');
    xlabel('Displacement (2 pi / 12)');
    ylabel('Correlation');
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
