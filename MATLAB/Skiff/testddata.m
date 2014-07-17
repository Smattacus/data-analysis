% testddata.m
% test out the mfile/function ddata.m
% 
% we can detect the power spectrum by direct fourier transform of the 
% time series.  Given two signals with the same coherent fluctuation
% but different photon statistics noise we can cross-correlate
% of course we can also simply add the light signals and this will 
% also lower the noise because of better statistics
%
% Note that the FFT of the xcorr differs from the square of the transform
% of the original signal by a factor of sqrt(NFFT).
%
NT=2^18;
t=1:NT;
t=t/NT;
freq=1500;
%y=3+.3*sin(2*pi*t*440);% perfectly coherent
% now consider a frequency that is not constant
y=3+.3*sin(2*pi*cumsum(freq*(1+randn(NT,1)))/NT);
%
R1=ddata(y,15);
R2=ddata(y,15);
R3=(R1+R2)/2;%normalized to have the same signal strength
[f h1]=spec(R1,1/NT);
[f g]=spec(y,1/NT);
[f h2]=spec(R2,1/NT);
[f h3]=spec(R3,1/NT);
figure(1); semilogy(f,abs(g).^2,f,abs(h1).^2,f,abs(h2).^2,f,abs(h3).^2)
xlabel('Frequency (Hz)')
title([num2str(freq) 'Hz tone' num2str(NT) 'sample freq average photon rate=3 amp=0.1'])
ylabel('transform.^2')
legend('sig','h1','h2','(h1+h2)/2')
%
y12=xcorr(R1-mean(R1),R2-mean(R2),'unbiased')*sqrt(NT);
y11=xcorr(R1-mean(R1),'unbiased')*sqrt(NT);
y22=xcorr(R2-mean(R2),'unbiased')*sqrt(NT);
y33=xcorr(R3-mean(R3),'unbiased')*sqrt(NT);
y12=y12(NT/2+(1:NT));
y11=y11(NT/2+(1:NT));
y22=y22(NT/2+(1:NT));
y33=y33(NT/2+(1:NT));
wind=((1-cos(2*pi*t))/2).^20';%this window has the same area under the curve as 1.0
%                    a filter with a width a few correlation lengths
%                    allows a kind of filtering from long time series
[f P11]=spec(y11.*wind,1/NT);
[f P22]=spec(y22.*wind,1/NT);
[f P33]=spec(y33.*wind,1/NT);
[f P12]=spec(y12.*wind,1/NT);
figure(2); semilogy(f,abs(P12),f,abs(P11),f,abs(P22),f,abs(P33))
xlabel('Frequency (Hz)')
title([num2str(freq) 'Hz tone' num2str(NT) 'sample freq average photon rate=3 amp=0.1'])
ylabel('transform.^2')
legend('P12','P11','P22','P(h1+h2)/2')
%%
%
err = sqrt(3);
ered = err * sqrt(sum(wind.^2) / (NT - 1));
figure(2); hold;
semilogy(f, ered + 0 * abs(P12), 'black');
%%
%Let's try propagating the errors for each one of Fred's correlation sets:
y12e = xcorr_err(R1 - mean(R1), R2 - mean(R2), sqrt(mean(R1)), sqrt(mean(R2)));
y11e = xcorr_err(R1 - mean(R1), R1 - mean(R1), sqrt(mean(R1)), sqrt(mean(R1)));
y22e = xcorr_err(R2 - mean(R2), R2 - mean(R2), sqrt(mean(R2)), sqrt(mean(R2)));
y33e = xcorr_err(R3 - mean(R3), R3 - mean(R3), sqrt(mean(R3)), sqrt(mean(R3)));
y12em = mean(y12e);
y11em = mean(y11e);
y22em = mean(y22e);
y33em = mean(y33e);
%%
%Remove outer half same way as earlier.
y12e = y12e(NT/2 + (1:NT));
y11e = y11e(NT/2 + (1:NT));
y22e = y22e(NT/2 + (1:NT));
y33e = y33e(NT/2 + (1:NT));
%%
figure(3); plot(t, y12e, t, y11e, t, y22e, t, y33e); xlabel('Time'); ylabel('Error bar value');
title('Error bars for different correlation arrays.');
legend('y12e', 'y11e', 'y22e', 'y33e');
%%
%With these means, let's see if that standard error prop works as expected:
y12em_s = sqrt(sum(wind.^2) / NT * y12em);
y11em_s = sqrt(sum(wind.^2) / NT * y11em);
y22em_s = sqrt(sum(wind.^2) / NT * y22em);
y33em_s = sqrt(sum(wind.^2) / NT * y33em);
figure(1); semilogy(f, abs(P12), 'b', f, abs(P12) * 0 + y12em_s, 'r'); xlabel('F(Hz)'); ylabel('log(|P|)');
title('P12 With error floor');
figure(2); semilogy(f, abs(P22), f, abs(P11) * 0 + y12em_s); xlabel('F(Hz)'); ylabel('log(|P|)');
title('P11 With error floor');
figure(3); semilogy(f, abs(P12), f, abs(P22) * 0 + y12em_s); xlabel('F(Hz)'); ylabel('log(|P|)');
title('P22 With error floor');
figure(4); semilogy(f, abs(P12), f, abs(P33) * 0 + y12em_s); xlabel('F(Hz)'); ylabel('log(|P|)');
title('P33 With error floor');