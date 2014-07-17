clear all; clc;

%% Initializing values

loc = 'D:\Lab Work\Data\7-5-2013\';                                         %Directory and file
list = dir(strcat(loc,'*.h5'));

Fs = 1E6;                                                              %Sampling freq
Fc = 1E5;                                                                   %Chopping freq

n = 5;                                                                     %number of files to analyze
chan = 32;                                                                   %number of channels to analyze
%% Open files

data = cell(n,1);                                                           %Initialize cell array for data

for i = 1:n
    data{i} = h5read(strcat(loc,list(i).name), '/PMT_DATA_8BIT');                              %Data is for data, H is for headers
end

%% Deleting Channels from unused PMT (If only 1 PMT was used)

PMT = 2;                                                                    % Select which PMT was used (1,2)

if PMT == 1
    for i = 1:n
    data{i}(chan/2+1:chan,:) = [];
    end
elseif PMT == 2
    for i = 1:n
    data{i}(1:chan/2,:) = [];
    end
end

%% Summing all channels
data_sum = cell(n,1);

for i = 1:n
    data_sum{i} = sum(data{i});                                                   %Sums all channels for each file and puts each file as a row in a matrix
end

%% Obtaining Phases
ph = cell(n,1);

for i = 1:n
    ph{i} = genPhase(data_sum{i},Fc,50,1/Fs);
end

clear data_sum;
%% Summing half and half channels
data_hsum = cell(n,1);

for i = 1:n
    data_hsum{i}(1,:) = sum(data{i}(1:chan/4,:));                                      %Sums all channels for each file and puts each file as a row in a matrix
    data_hsum{i}(2,:) = sum(data{i}(chan/4+1:chan/2,:));
end

clear A;
%% Generating Top and Bottoms
top = cell(n,1);bot = cell(n,1);

for i = 1:n
    [top{i}(1,:), bot{i}(1,:)] = getTopBot2(data_hsum{i}(1,:),...
        square(pi/2+ph{i}),Fs,Fc/2);
    [top{i}(2,:), bot{i}(2,:)] = getTopBot2(data_hsum{i}(2,:),...
        square(pi/2+ph{i}),Fs,Fc/2);
end

%% Substracting Bottoms from Tops
diff_tb = cell(n,1);

for i = 1:n
    diff_tb{i}(1,:) = top{i}(1,:)-bot{i}(1,:);
    diff_tb{i}(2,:) = top{i}(2,:)-bot{i}(2,:);
end
    
%% Cross Correlations of half sums
Axh = cell(n,1);

for i = 1:n
    Axh{i} = xcorr(diff_tb{i}(1,:)-mean(diff_tb{i}(1,:)),...                    %Cross correlation between files, using first file as reference
        diff_tb{i}(2,:)-mean(diff_tb{i}(2,:)), 'unbiased');
end

%% Averaging Cross correlations
c = length(Axh{1});
Axh_avg = zeros(1,c);

for i = 1:n
    Axh_avg = Axh_avg+Axh_top{i};
end

Axh_avg = Axh_avg/n;

%% Checking correlation length

N = (c+1)/2;
tt = (-(N-1):(N-1))/Fc;
figure(1); plot(tt,Axh_avg);

%% Windowing

tc = .1;                                                                    %correlation time
t = linspace(-1,1,c);                                                       %Create Window
win = exp(-t.^2/(tc^2/2));

Axh_win = Axh_avg.*win;                                                     %Multiply cross correlations by window

%% Fourier Transform
[f, g] = spec(Axh_win, 1/Fc);

%% Plot
figure(1); plot(f,abs(g));
xlabel('Frequency (Hz)'); ylabel('Power');
title('Demodulated Signal (half-half cross corr)');