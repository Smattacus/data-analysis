%%
%Script to generate arrays of average LIF photons for the Zeeman data from
%11/5/2013.
list = ls('*Zeeman*.h5');
Nf = size(list,2);
[T1, B1, T2, B2] = genAllTPBoxFindPhase(p, 1e5, 1e6);
avg_photons = mean(T1 - B1);