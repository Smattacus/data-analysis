%Script to print phases to a file in MATLAB.
%Template
xcmpath = '/Users/smattingly/Data/11-5-2013/XCMeans';
dyename = 'dyehere';
datapath = sprintf('/Users/smattingly/Data/11-5-2013/V2PT/%s', dyename);
savename = sprintf('%s/%s.mat', xcmpath, strrep(dyename, '.', '_'));

saveXCMeanBoxFindPhase(datapath, savepath, 1e5, 1e6);