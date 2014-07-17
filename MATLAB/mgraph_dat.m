function [] = mgraph_dat(data)
A = data;
AA = double(A);
sA = size(AA);
dim = sA(1) / 32;
MA = reshape(AA, 32, dim);
mesh(1:dim, 1:32, MA);
n = strrep(filename, '_', ' ')
title(['Mesh plot for ', n]);
xlabel('Time');
ylabel('PMT #');