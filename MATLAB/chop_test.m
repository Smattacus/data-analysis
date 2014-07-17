%%Test dataset for demodulation of a signal
Fs = 100000;
t = 0:.00001:0.99999;
A = 5 + sin(2 * pi * 180 * t) + .5 * cos(2 * pi * 60 * t) + rand(1,Fs);
B = rand(1, Fs);
C = chop_modulate(A, B, 50);
D = chop_modulate(A, B, 20);
%C is the modulated signal - it is "on" and "off" every 50 elements.
%e.g. 20 Hz signal chopping.
%Perform fourier analysis:
Y = fft(C);
L = size(C);
x = Fs/2 * linspace(0, 1, L(2) / 2);
Yd = fft(D);
Ld = size(D);
xd = Fs/2 * linspace(0, 1, L(2) / 2);

%Create a phase - varying signal
temp = zeros(1, 150000);
for i = 1:max(size(t)) / 1500;
    a1 = splice(A(1 + (i-1) * 1500:(i-1) * 1500 + 500), ...
        B(1 + (i-1) * 1500:(i-1) * 1500 + 500), 50);
    a2 = splice(A(501 + (i-1) * 1500: (i-1) * 1500 + 1020), ...
        B(501 + (i-1) * 1500: (i-1) * 1500 + 1020), 52);
    a3 = splice(A(1021 + (i-1) * 1500:(i-1) * 1500 + 1500), ...
        B(1021 + (i-1) * 1500:(i-1) * 1500 + 1500), 48);
    temp(1 + (i-1) * 1500:i*1500) = [a1 a2 a3];
end
E = temp;

%Create a different phase varying signal.
Fs2 = 100000;
t2 = 0:.00001:.99999;
L2 = max(size(t));
A2 = 5 + sin(2 * pi * 180 * t) + 1.5 * cos(2 * pi * 60 * t) + 1/2 *rand(1,L2);
B2 = 1/2 * rand(1, L2);
C2 = splice(A2(1:L2/2), B2(1:L2/2), 50);
D2 = splice(A2(L2/2 + 1: L2), B2(L2/2 + 1:L2), 51);
F = [C2 D2];