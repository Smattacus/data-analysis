function bitArr = breakPmt8Bit(A)
%
%bitArr = breakPmt8Bit(A)
%
%Function to take an 32xN of 4byte characters and decompose it
%into a 32x4N array of 1byte PMT readings.
%
%This functions ASSUMES that the STRUCK SIS3820 is set to
%acquire data in 8bit mode.
%
%

%Construct mask array for bits 1 & 2
mask = uint32(255); %equivalent to 0x000000FF
%Construct final array
bitArr = uint8(zeros(size(A,1), size(A,2) * 4));
ba_c = size(bitArr,2);
sa = size(A,1) * size(A,2);
bitArr(1:4:end) = uint8(bitand(A, mask));
bitArr(2:4:end) = uint8(bitand(bitshift(A, -8), mask));
bitArr(3:4:end) = uint8(bitand(bitshift(A, -16), mask));
bitArr(4:4:end) = uint8(bitand(bitshift(A, -24), mask));

end
%{
for i=1:4:(size(A, 1)* size(A,2))
    for k=1:4
        e = [];
        for j=1:4
            t = bitshift(A(i + k), - (j * 8));
            t = bitand(t, mask);
            e = [e; e];
        end
        bitArr = [bitArr, e];
    end
end
%}