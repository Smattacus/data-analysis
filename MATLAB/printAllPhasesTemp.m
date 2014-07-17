%Script to print phases to a file in MATLAB.
%Template
path = '/Users/smattingly/Data/11-5-2013/V2PT/Dye_611.65975'
cd path
mkdir Phase_PNG

list = dir('*h5');
Nf = size(list,1);
for i=1:Nf
    fn = strrep(list(1).name, 'num0', sprintf('num%d', i));
    findMaxPhaseFilePlot(fn, 8, 1e6, 1e5, true);
end