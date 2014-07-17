function bitArr = breakPmtBits(A)
%Function to take an array of 4byte characters and decompose it
%into a 32xN array of 1byte PMT readings.
%
%Bitmasks to use. PMTs are ordered from right to left.
%
%Construct mask array
mask = 3 * ones(size(A,1), 1);
%Construct final array
bitArr = zeros(size(A,1)/8, 32);
A = bitand(A())