%
%Script to generate a 2D surface plot of spectra vs. wavelength (i.e.
%velocity)

%%
%Generate the velocities.
importDyes('dye_areas.txt');

dyes = data(:,1);
dye_lambda0 = 611.6598968649599;
diode_lambda0 = 668.6113026856146;
diode = 668.61076;
v_diode = dopplerVelocity(diode_lambda0, diode);
v_dyes = dopplerVelocity(dye_lambda0, dyes);

%%
%Get the areas from the data
areas1 = data(:,3);
areas2 = data(:,5);

%Plot the un - normalized correlation function
dye_axis = dyes - 611.6;
figure(1); plot(v_diode - v_dyes, areas1); xlabel('Diode - Dye Velocity (m/s)'); ylabel('Area under fitted Lorentzian'); title('Small Peaks');
figure(2); plot(v_diode - v_dyes, areas2); xlabel('Diode - Dye Velocity (m/s)'); ylabel('Area under fitted Lorentzian'); title('Large Peaks');
figure(3); plot(v_diode - v_dyes, areas1 + areas2); xlabel('Diode - Dye Velocity (m/s)'); ylabel('Area under fitted Lorentzian'); title('Sum of Areas');

%%
%
list = dir('xcmean_*.mat');
N = 800000;
t = (-(N-1):(N-1))/1e5;
tc = 0.01;
win = exp(-(t/tc).^2/2);

load(list(1).name);
[f, g] = spec(xcmean, 1e-5);
ifg = find(abs(f) < 5000);
size(ifg)
specs = zeros(size(ifg,2), size(list,1)-1);

for i=1:(size(list,1)-1)
    load(list(i).name);
    [f, g] = spec(xcmean .* win, 1e-5);
    %Try subtracting the means:
    ig = find((f < 3e4) - (f < 3e3));
%    g = abs(g) - mean(abs(g(ig)));
    g = abs(g);
    specs(:,i) = g(ifg);
end

mesh(dye_axis, f(ifg), specs); xlabel('Dye Wavelength (nm + 611.6)');
ylabel('Frequency'); colorbar; title('Spectra vs. Wavelength Surface Plot');
view(0, 90)

%%
%Try again with a sorted list - make sure it looks the same.
newlist = list(newI);
for i=1:(size(newlist,1)-1)
    load(newlist(i).name);
    [f, g] = spec(xcmean .* win, 1e-5);
    specs(:,i) = g(ifg);
end

mesh(f(ifg),dye_axis, abs(specs).'); ylabel('Acquisition Order (num)');
xlabel('Frequency'); colorbar; title('Spectra vs. Wavelength Surface Plot');
view(0, 90)

%%
%Plot the Lines as a function of time.
importfile_onlynums('time_list_dechours.txt');
dye_times = data(:,1);
hours = data(:,2);
[Sdyes, I] = sort(dye_times);
unsorted = 1:length(Sdyes);
newI(I) = unsorted;
h = waterfall(f(ifg), hours, abs(specs(:, newI).')); view(0, 90); xlabel('Frequency (Hz)');
ylabel('Acquisition Time (hours)'); title('Waterfall Plot of Spectra vs. Time');
set(h, 'LineWidth', 1.5);
