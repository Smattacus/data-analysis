function out=ddata(na,nm)
%
%function to create simulated fluctuation data given the time series na
%
%     out=ddata(na,nm)
%
%na is the instantaneous average photon flux which determines the
%probability for photon counting at each time step (Poisson statistics)
%there is probably a clever way to determine how high the count numbers can
%go (nm) but this is left as an input parameter here.  A saturation of
%out at nm means that you have not set nm high enough.
%out is a time series with the same length as na and statistically with the
%same magnitude.
%
LY=length(na);
h=rand(LY,1);%Output has the same size as input array.
LN=nm;%someday find clever way to determine the best size of nm
[N,NA]=meshgrid(0:(LN-1),na);
G=[zeros(LY,1) cumsum(exp(-NA).*(NA.^N)./factorial(N),2)];
%rand gives you an equal probability between 0->1 so I integrate
%the Poisson probability vs n using cumsum and see where this integral 
%exceeds the output of rand.  The index of this point is the random number
%of photons for the time stip.  I use diff to single out the first point
%where the integral exceeds the random number and then a vector
%multiplication to find the index of that "zero crossing".
%this result is the simulated number of photons for each time step.
T=.5*diff(sign(G-h*ones(1,LN+1)),[],2);%Identifies the first nonzero.
NN=(0:(LN-1))';
out=T*NN;