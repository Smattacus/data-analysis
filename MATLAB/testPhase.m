clear all; clc;

%% Initializing values
%Directory and file
loc = '/Users/Sean/Documents/Skiff/Data/11-5-2013/LIFTest/';                                         
list = dir(strcat(loc,'*Garray*.h5'));
%list = dir('*GAGE_LIF*.h5');
%Sampling freq
Fs = 1E6;  
%Chopping freq
Fc = 1E5;                                                                   
%Number of files to analyze
n = 3;                     
%PMT used (1 = 1, 2 = 2, 3 = both)
PMT = 3;
%Acquisition Time
Ta = 8;

%% Open files
%Data for h5 files
data = h5read(strcat(loc,list(1).name), '/PMT_DATA_8BIT');

%% Generating Phase for each PMT
data_sum(1,:) = sum(data(7:9,:));
data_sum(2,:) = sum(data(23:25,:));
%Phase generated from code
ph(1,:)= genPhase(data_sum(1,:),Fc,50,1/Fs);
ph(2,:)= genPhase(data_sum(2,:),Fc,50,1/Fs);
%Regular phase for square wave with freq Fc
ph2 = linspace(0,2*pi*Fc*Ta,Fs*Ta);
p = linspace(0,2*pi,100);

%% Obtaining LIF from Tops and
LIF = zeros(2,100); LIF2 = zeros(2,100);
%top = zeros(2,8E6); top2 = top; bot = top;

for i = 1:100
    top(1,:) = data_sum(1,:).*(square(p(i)+ph(1,:))+1)/2;
    top(2,:) = data_sum(2,:).*(square(p(i)+ph(2,:))+1)/2;
    top2(1,:) = data_sum(1,:).*(square(p(i)+ph2)+1)/2;
    top2(2,:) = data_sum(2,:).*(square(p(i)+ph2)+1)/2;
    LIF(1,i) = sum(top(1,:));
    LIF(2,i) = sum(top(2,:));
    LIF2(1,i) = sum(top2(1,:));
    LIF2(2,i) = sum(top2(2,:));
    disp(i);
end

%% Generating plots
figure(3); plot(p,LIF(1,:)); title('PMT1 Function Generated Phase 0%AM');
figure(4); plot(p,LIF(2,:)); title('PMT2 Function Generated Phase 0%AM');
figure(5); plot(p,LIF2(1,:)); title('PMT1 Regular Phase 0%AM');
figure(6); plot(p,LIF2(2,:)); title('PMT2 Regular Phase 0%AM');