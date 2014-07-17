%Test for cell mode structure generating plots
%%
t = (1:2000)/(2 * pi * 50);
for i=1:4
%%
    s = sin(t + i * pi / 4); 
    plot(t, s);
end