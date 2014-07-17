%%
%Script for experimenting with demodulating
%This is demodulation using FFT and IFFT operations,
%NOT Hilbert transforms.
A = 1;
tp = 2^12;
omega1 = 240;
omega2 = 10;
gamma = pi/7;
beta = 5;
t = 0:2*pi / tp:2*pi*(1-1/tp);
x = A * cos(omega1*t + gamma + beta*sin(omega2*t));
plot(t(1:round(length(x)/omega2)), x(1:round(length(x)/omega2)))
xlim([0 t(round(length(x)/omega2))]);
xlabel('time (s)')
ylabel('magnitude')
%%
%FFT
ff = fft(x);
ax = linspace(-tp/pi/4, (tp-2)/pi/4,tp);
%db display.
dbf = 20 * log10((abs(ff) / length(ff)*2 + 10e-12)/10e-12);
dbf = fftshift(dbf);
%Manipulate axes
plot(ax, dbf);
xlim([-tp/pi/4 tp/pi/4]);
ylim([0 1.2*max(dbf)])
xlabel('frequency (rad/s)')
ylabel('magnitude dB rel. e-12');
%Set negative freqs to zero, double positive freqs.
%%
%Delete negative freqs, double positive freqs. Leave 0.
gf = ff;
gf(2:end) = 2 * ff(2:end);
gf(end/2+1:end) = 0;
g = ifft(gf);
dbg = 20 * log10((abs(gf)/length(gf)*2 + 10e-12)/10e-12);
dbg = fftshift(dbg);
plot(ax, dbg);
xlim([-tp/pi/4 tp/pi/4]);
ylim([0 1.2 * max(dbg)]);
xlabel('Frequency (rad/s)');
ylabel('Magnitude dB rel. 10^-12');

%%
%g is now the analytic signal. Instead of doing all this,
%we could have said g = hilbert(x).
%Find the phase info:
pha = angle(g);
upha = unwrap(pha);
%%
%This is to take care of Omega * t + gamma
p = polyfit(t, upha, 1);
p(2) = upha(1);
omega1t = polyval(p,t);
uphas = upha - omega1t; %linear offset sub
%FFT the other stuff:
mf = fft(uphas); %Spectrum of phase demo signal
dbmf = 20 * log10((abs(mf)/length(mf)*2 + 10e-12)/10e-12);
dbmf = fftshift(dbmf);
plot(ax, dbmf);
xlim([0 10 * omega2/2/pi]);
ylim([min(dbmf) 1.2*max(dbmf)]);