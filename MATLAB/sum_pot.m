function [ num ] = sum_pot( n, m )
%SUM Summary of this function goes here
%   Detailed explanation goes here

num = 32 / pi^2 * 1 / (n * m) * sinh(pi / 2 * sqrt(n^2 + m^2)) / sinh(pi * sqrt(n^2 + m^2)) ;
end

