clear all
close all
clc

load('Data_Rx');

frame_len = length(Data_10dB(:,1,1));
numb_symb_len = length(Data_10dB(1,1,:));
symb_len = 80;

Data_1dB = reshape(Data_1dB,frame_len*numb_symb_len,1);
Data_5dB = reshape(Data_5dB,frame_len*numb_symb_len,1);
Data_10dB = reshape(Data_10dB,frame_len*numb_symb_len,1);
Data_15dB = reshape(Data_15dB,frame_len*numb_symb_len,1);
Data_20dB = reshape(Data_20dB,frame_len*numb_symb_len,1);

mean_1dB = mean(Data_1dB);
mean_5dB = mean(Data_5dB);
mean_10dB = mean(Data_10dB);
mean_15dB = mean(Data_15dB);
mean_20dB = mean(Data_20dB);

std_1dB = std(Data_1dB);
std_5dB = std(Data_5dB);
std_10dB = std(Data_10dB);
std_15dB = std(Data_15dB);
std_20dB = std(Data_20dB);

var_1dB = var(Data_1dB);
var_5dB = var(Data_5dB);
var_10dB = var(Data_10dB);
var_15dB = var(Data_15dB);
var_20dB = var(Data_20dB);

%% graficas
subplot(1,5,1);
x = histogram(real(Data_1dB));
xlabel("Mean = " + mean_1dB)
text(mean_1dB,max(x.Values),'x Mean')

subplot(1,5,2);
x = histogram(real(Data_5dB));
xlabel("Mean = " + mean_5dB)
text(mean_5dB,max(x.Values),'x Mean')

subplot(1,5,3);
x = histogram(real(Data_10dB));
xlabel("Mean = " + mean_10dB)
text(mean_10dB,max(x.Values),'x Mean')

subplot(1,5,4);
x = histogram(real(Data_15dB));
xlabel("Mean = " + mean_15dB)
text(mean_15dB,max(x.Values),'x Mean')

subplot(1,5,5);
x = histogram(real(Data_20dB));
xlabel("Mean = " + mean_20dB)
text(mean_20dB,max(x.Values),'x Mean')

%% PEAK TO AVARAGE

symb = Data_1dB(1:symb_len);
acc = 0;
max_symb = max(abs(symb));
mean_symb = mean(abs(symb));

peak_to_avarage = (max_symb^2)/(mean_symb^2);
ptoa_db = mag2db(peak_to_avarage);
