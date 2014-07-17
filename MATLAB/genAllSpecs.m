%Script for generating all of the PMT Spectrograms in order to determine
%good windows for the chop frequency (or if a file should be rejected since
%the LIF signal is too weak.
%%
list = dir('*.DAT');
for i=1:size
%%
    fn = list(i).name;
    plotPMTSpec(fn);
end