function P=poisson(N,na)
%poisson
%Calculate the poisson distribution P(N,na) given
%   N the Maximun photon number  (return P for n=0,1, ... N)  N natural
%   number.
%   na the average photon number, real, positive number
%   P=poisson(N,na)
%   na can be a 1-D vector
%   P has size (N+1) x length(na)
%
[nna,n]=meshgrid(na(1:end),1:N);
[na0,n0]=meshgrid(na(1:end),0);
P=[exp(-na0); exp(-cumsum(log(n))-nna +n.*log(nna))];