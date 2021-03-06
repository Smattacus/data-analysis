t=(1:44100)/44100;
y=3+.1*sin(2*pi*440*t);
LY=length(y);
h=rand(LY,1);%Output has the same size as input array.
LN=10;%find clever way to determine the best size of N
[N,NA]=meshgrid(0:(LN-1),y);
G=[zeros(LY,1) cumsum(exp(-NA).*(NA.^N)./factorial(N),2)];
T=.5*diff(sign(G-h*ones(1,LN+1)),[],2);%Identifies the first nonzero.
NN=(0:(LN-1))';
out=T*NN;