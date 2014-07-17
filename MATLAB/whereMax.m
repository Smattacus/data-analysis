function freq_max = whereMax(f, g)
%Quick function to determine the location of the maximum near 100 kHz
%
%   freq_max = whereMax(f, g)
%
ig = find((f < 105000) - (f < 95000));
temp = f(ig);
i = find(max(abs(g(ig))) == abs(g(ig)));
freq_max = temp(i);