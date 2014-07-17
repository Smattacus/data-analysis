function [f, mat] = mesh_onechan(directory, channum)
%This is a function for displaying and returning the mesh array for all of
%the graphs for a particular channel. It assumes that the arrays have been
%saved with the filename convention Dye_*.mat, and that each file has the
%ac_ch1_mn, ac_ch1_std, ac_ch2_mn, ac_ch2_std variables.
%
%[f, mat] = mesh_onechan(directory, channum)
%
%INPUTS:
%   directory   - Directory the files are in.
%   channum     - 1 or 2. Any other entry will result in error return.
%OUTPUTS:
%   f           - Frequency of transform
%   mat         - Matrix of transforms.
%
%This function was made for 5-30-2013 and later data.
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
   if channum == 1
       [f, gw] = spec((ac_ch1_mn - mean(ac_ch1_mn)) .* win, 1e-5);
   elseif channum == 2
       [f, gw] = spec((ac_ch2_mn - mean(ac_ch2_mn)) .* win, 1e-5);
   else
       return
   end
   mat(:,i) = gw;
end
ig = find((f < 2500) - (f < 1000));
mesh(-(1:Nfiles), f(ig), abs(mat(ig, :))); view(0,90); xlabel('Num');
ylabel('Frequency'); title(sprintf('Surface plot of Data from channel %d, \n%s', channum, directory));