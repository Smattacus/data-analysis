function [C] = vary_phase(A, B, n);
%Splices two vectors on top of each other BUT it varies the phase
%throughout the signal.
%Assume the vectors are the same size.
%Inputs:
%   A - First vector
%   B - Second vector
%   n - Number of elements
%The phase is varied as follows: every 10 periods, the phase is adjusted.
%It starts with a period of n, then n+2, then n-2. These seemingly
%arbitrary numbers are to test the hilbert transform to recreate
%the phase variation of the chop.
 L = max(size(A));
 temp = zeros(1, L);
 for i:
     
 end
 
 function [C] = splice(A, B, n);
%Splices two vectors on top of each other in order to simulate
%chopper modulation.
%Assume the vectors are the same size.
%Inputs:
%   A - First vector
%   B - Second vector
%   n - number of elements to copy from B into A.
%The output vector will be of the form
% C = [A(1:n) B(n+1:2n) A(2n+1:3n) ...];
%Assume they are column vectors
L = size(A);
temp = zeros(1, L(2));
for i=1:L(2)/(2*n)
    temp(1 + (i-1) * 2*n:2*n*i) = [ones(1,n) zeros(1,n)];
end
C = (temp == 1) .* A + (temp == 0) .* B;
 end