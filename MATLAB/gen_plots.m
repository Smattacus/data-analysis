%Script to generate plots from the file
%plot_suffix_list.txt

f = fopen('plot_suffix_list.txt', 'r');
format = '%s %s';
c = textscan(f, format);

%Column array of strings of paths and suffixes
paths = c{1};
suffixes = c{2};

for i=0:size(suffixes,1)
   p = strtrim(paths(i));
   suff = strtrim(suffixes(i));
   suff_strip = strrep(suff, '.txt');
   fn_n = sprintf('Nplot_run_%s.jpg', suff_strip);
   fn_phi = sprintf('PhiPlot_run_%s.jpg', suff_strip);
   fn_T = sprintf('TPhiPlot_run_%s.jpg', suff_strip);
   plotRadPhi(p, suff);
   print(fn_n, '-djpeg');
   plotRadT(p, suff);
   print(fn_phi, '-djpeg');
   plotRadialDensity(p, suff);
   print(fn_T, '-djpeg');
end