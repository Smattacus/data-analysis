%FFT play script
fs = 1000;
t = 0:.001:.999;
y = 1.5 * sin(2 * pi * 50 * t) + sin( 2 * pi * 120 * t) + rand(1, 1000);
plot(fs * t(1:100), y(1:100));
xlabel('Time in ms');
ylabel('Sampling');
Y = fft(y);
s = size(Y);
f = fs / 2 * linspace(0, 1, s(2)/2);