function pb = genBasePhase(total_t)
%FUNCTION pb = genBasePhase(total_t)
%
%This function generates a base phase array
%for a <total_t> second long file which is oversampled at 10x.
pbone = linspace(0, 2 * pi - 2 * pi / 10, 10);
pb = [];
%This is stupid and ugly, but faster than one loop by far.
for i=1:1000
    pb = [pb, pbone];
end
pbone = pb;
pb = [];
for i=1:50
    pb = [pb, pbone];
end
pbone = pb;
pb = [];
for i=1:total_t*2
    pb = [pb, pbone];
end