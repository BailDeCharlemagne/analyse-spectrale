clear all; 
close all;
pkg load control;
pkg load signal; 
addpath(genpath('/home/mony/octave/'));
 


%% Transform flute signal to data

[flute, fs] = audioread('fluteircam.wav', native_float_format);
flute_t = transpose(flute);
t_ech = 1:length(flute);
t = t_ech/fs; %échelle temporelle
f_subsampling = 4;
figure;
plot(t, flute);
title('Signal pur de la flute - non bruité'); 

%%%% Add Noise to Flute
%%%
%%%Eb_No = 0.05; %en dB 
%%%Eb = (norm(flute)^2)/length(flute);
%%%No = Eb_No*log10(Eb); 
%%%variance = No/2; %Revoir avec Telecom
%%%bruit = randn(length(flute),1)*variance; 
%%%noisy_flute = bruit + flute;
%%%figure; 
%%%plot(t, noisy_flute) 
%%%title('Signal bruité de la flute');  
%%%
%%%
%%%%% FFT du signal 
%%%
%%%figure;
%%%hold on; 
%%%
%%%noisy_fft_complex = abs(fft(noisy_flute)); 
%%%noisy_fft_real = noisy_fft_complex(1:length(noisy_flute)/2); 
%%%freq_red = fs*(0:length(noisy_flute)/2-1)/length(noisy_flute);
%%%plot(freq_red, noisy_fft_real, 'r');
%%%
%%%y_fft_complex = abs(fft(flute_t)); 
%%%y_fft_real = y_fft_complex(1:length(flute)/2); 
%%%freq_red = fs*(0:length(flute)/2-1)/length(flute);
%%%plot(freq_red, y_fft_real, 'b');
%%%hold off;
%%%
%%%audiowrite('fluteircam-bruité.wav', noisy_flute, fs);
%
%%%%%% On s'intéresse à la séparation haute fréquence entre 4 et 6 secondes.
%
%close all; 
%t_debut = 4; %en sec 
%t_fin = 6; %en sec
%t_debut_ech = t_debut*length(flute)/max(t); 
%t_fin_ech = t_fin*length(flute)/max(t);  
%flute_morceau = flute_t(t_debut_ech:t_fin_ech); 
%
%% Pisarenko
% 
%[ff_pisa, mydsp_pisa] = mypisarenko(flute_morceau, 100, fs, 1);
%title('Pisarenko'); 
%% On ne veut prendre que la partie où les fréquences sont positives 
%ff_pisa_real = ff_pisa(floor(length(ff_pisa)/2)+1:end); 
%mydsp_pisa_real = mydsp_pisa(floor(length(mydsp_pisa)/2)+1:end);
%plot(ff_pisa_real, mydsp_pisa_real);
%[pks_pisa loc_pisa] = findpeaks(mydsp_pisa_real); 
%index_pks_pisa = ff_pisa_real(loc_pisa)

% test échantillonnage 

