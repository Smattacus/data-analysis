function [f, mat] = mesh_xcmeans(directory)
%Function to generate the mesh plot of a directory of xcmean files.
%Windowing is used, with 0.01 tc.
%
%mesh_xcmeans(directory)
%
list = dir(sprintf('%s/*-means.mat', directory));
Nfiles = size(list,1);
load(list(1).name);
N = size(xcmean,2);
N = (N+1)/2;
t = (-(N-1):(N-1))/1e5;
win = exp(-(t/.01).^2/2);
mat = zeros(2 * N - 1, size(list,1));
for i=1:Nfiles
   load(list(i).name);
   [f, gw] = spec((xcmean - mean(xcmean)) .* win, 1e-5);
   mat(:,i) = gw;
end
ig = find((f < 2500) - (f < 1000));
mesh(-(1:Nfiles), f(ig), abs(mat(ig, :))); view(0,90); xlabel('Num');
ylabel('Frequency'); title(sprintf('Surface plot of Data from \n%s', directory));