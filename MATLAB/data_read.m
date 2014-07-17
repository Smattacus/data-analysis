%Read the first second.
f1 = 'plasma_noise_lights_Tue_Apr_17_23_18_57_pmt_800mV_s0_num0.dat';
f2 = 'plasma_noise_lights_Tue_Apr_17_23_18_57_pmt_800mV_s0_num1.dat';
f3 = 'plasma_noise_lights_Tue_Apr_17_23_18_57_pmt_800mV_s0_num2.dat';
f4 = 'plasma_noise_lights_Tue_Apr_17_23_18_57_pmt_800mV_s0_num3.dat';
f5 = 'plasma_noise_lights_Tue_Apr_17_23_18_57_pmt_800mV_s0_num4.dat';

%Readout all the data
[A1,HEAD] = hread_PMTdata(f1);
[A2] = read_PMTdata(f2);
[A3] = read_PMTdata(f3);
[A4] = read_PMTdata(f4);
[A5] = read_PMTdata(f5);

sizes = [size(A1); size(A2); size(A3); size(A4); size(A5)];

data1 = reshape(A1, 32, sizes(1,1) / 32);
%data2 = reshape(A2, 32, sizes(2,1) / 32);
%data3 = reshape(A3, 32, sizes(3,1) / 32);
%data4 = reshape(A4, 32, sizes(4,1) / 32);
%data5 = reshape(A5, 32, sizes(5,1) / 32);

pmt1 = double(data1(1,:));
pmt4 = double(data1(4,:));
pmt6 = double(data1(6,:));

pmt1 = pmt1 - mean(pmt1);
pmt4 = pmt4 - mean(pmt4);
pmt6 = pmt6 - mean(pmt6);

ac1 = xcorr(pmt1);
ac4 = xcorr(pmt4);
ac6 = xcorr(pmt6);

%use the psd / spectrum stuff:
%Hamming window, periodic sampling, 385 element length (which
% corresponds to ~ 2 milliseconds.)
h = spectrum.welch({'hamming', 'symmetric'}, 385);
psd(h, ac1, 'Fs', 192307.692);
hold
psd(h, ac4, 'Fs', 192307.692);
psd(h, ac6, 'Fs', 192306.692);