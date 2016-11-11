clear all; 
close all;
pkg load control;
pkg load signal; 
addpath(genpath('/home/mony/octave/'));
 

%%%% Transform flute signal to data

[flute, fs] = audioread('fluteircam.wav', native_float_format);
flute_t = transpose(flute);
t_ech = 1:length(flute);
t = t_ech/fs; %échelle temporelle
f_subsampling = 4;
%figure;
%plot(t, flute);
%title('Signal pur de la flute - non bruité'); 



% Add Noise to Flute

Eb_No = 2; %en dB 
Eb = abs((norm(flute)^2))/length(flute);
No = Eb_No*log(-Eb/10); 
variance = No/2; %Revoir avec Telecom
bruit = randn(length(flute),1)*variance; 
noisy_flute = bruit + flute;
figure; 
plot(t, noisy_flute) 
title('Signal bruité de la flute');  


%% FFT du signal 

figure;
hold on; 

noisy_fft_complex = abs(fft(noisy_flute)); 
noisy_fft_real = noisy_fft_complex(1:length(noisy_flute)/2); 
freq_red = fs*(0:length(noisy_flute)/2-1)/length(noisy_flute);
plot(freq_red, noisy_fft_real, 'r');

y_fft_complex = abs(fft(flute_t)); 
y_fft_real = y_fft_complex(1:length(flute)/2); 
freq_red = fs*(0:length(flute)/2-1)/length(flute);
plot(freq_red, y_fft_real, 'b');
hold off;

audiowrite('fluteircam-bruité.wav', noisy_flute, fs);
