function [GcosM, GsinM] = createGMatrix(varargin)
%
% G = createGMatrix(path, suffix)
% G = createGMatrix()
%
%No Inputs: Prompts the user for a file.
%2 Inputs:
% path      - Path to the desired file.
% suffix    - Suffix on the desired data run.
%Outputs:
% G         - G matrix ready to be used in an SVD algorithm

if (nargin == 2) 
    list = getArrayFileList(varargin{1}, varargin{2});
else
    list = getArrayFileList();
end

size(list)

GsinM = []; GcosM = [];
for i=1:size(list,1)
   filename = list(i).name;
   [angle, Gsin, Gcos] = genSVDComp(filename);
   size(Gsin);
   size(Gcos);
   GsinM = [GsinM, Gsin];
   GcosM = [GcosM, Gcos];
end
