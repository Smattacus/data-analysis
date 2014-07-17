%%Cell mode script to generate a normalized Zeeman effect plot.
%%
load files
N = size(lefts, 1);
dyes = [];
diodes = [];
%%
for i=N
    [A, H] = getPMTData(fns{i}, true);
    

%%
for i=1:N
    display(sprintf('Working on %s', fns{i}))
    [A, H] = getPMTData(fns{i}, true);
    [diode_n, dye_n] = get_LIFNormStrength(A, lefts(i), rights(i), 0.01);
    dyes = [dyes, dye_n];
    diodes = [diodes, diode_n];
end

