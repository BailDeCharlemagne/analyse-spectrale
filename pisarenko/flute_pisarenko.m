%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REGLAGES PRELIMINAIRES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;
pkg load control;
pkg load signal; 
addpath(genpath('/home/mony/octave/'));
 

%%% Transform flute signal to data

[flute, fs] = audioread('fluteircam.wav', native_float_format);
flute_t = transpose(flute);
t_ech = 1:length(flute);
t = t_ech/fs; %échelle temporelle
f_subsampling = 4;
%figure;
%plot(t, flute);
%title('Signal pur de la flute - non bruité'); 

%%% Add noise to signal 

Eb_No = 2; %en dB 
Eb = abs((norm(flute)^2))/length(flute);
No = Eb_No*log(-Eb/10); 
variance = No/2; %Revoir avec Telecom
bruit = randn(length(flute),1)*variance; 
noisy_flute = bruit + flute;
figure; 
%plot(t, noisy_flute) 
%title('Signal bruité de la flute');  

%%% Split flute signal into 30ms intervals 

interval = 30*10^(-3); 
interval_ech = interval*fs;
nbIntervals = floor(length(flute_t)/interval_ech); 
flute_ech = zeros(nbIntervals, interval_ech);
for i=1:nbIntervals
  flute_ech(i,:) = flute_t((i-1)*interval_ech+1:i*interval_ech); 
end; 


###########################################################################
##########            PISARENKO SUR 4 A 6 SECONDES              ###########
###########################################################################

%%% On s'intéresse à la séparation haute fréquence entre 4 et 6 secondes.

close all; 
t_debut = 4; %en sec 
t_fin = 6; %en sec
t_debut_ech = t_debut*length(flute)/max(t); 
t_fin_ech = t_fin*length(flute)/max(t);  
flute_morceau = flute_t(t_debut_ech:t_fin_ech); 

%%% Pisarenko sur le signal pur 

# Pisarenko 0
[ff_pisa_0, mydsp_pisa_0] = mypisarenko(flute_morceau, 100, fs, 0); % donne densité de puissance spectrale
ff_pisa_real_0 = ff_pisa_0(floor(length(ff_pisa_0)/2)+1:end); 
mydsp_pisa_real_0 = mydsp_pisa_0(floor(length(mydsp_pisa_0)/2)+1:end);
figure; 
plot(ff_pisa_real_0, mydsp_pisa_real_0);
title('Paramètre 0') 

# Pisarenko 1 
[ff_pisa1, mydsp_pisa1] = pisarenko_real(flute_morceau, 100, fs, 1, 1); % donne spectre de raies
% On cherche un vecteur de pics de fréquences 
[pks_pisa1 loc_pisa1] = findpeaks(mydsp_pisa1); 
index_pks_pisa1 = ff_pisa1(loc_pisa1)

# Pisarenko 2 : 
[ff_pisa2, mydsp_pisa2] = pisarenko_real(flute_morceau, 100, fs, 1,1); % donne enveloppe interpolée
% On cherche un vecteur de pics de fréquences 
[pks_pisa2 loc_pisa2] = findpeaks(mydsp_pisa2); 
index_pks_pisa_2 = ff_pisa2(loc_pisa2) 
 

