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

%%%%% On s'intéresse à la séparation haute fréquence entre 4 et 6 secondes.

close all; 
t_debut = 4; %en sec 
t_fin = 6; %en sec
t_debut_ech = t_debut*length(flute)/max(t); 
t_fin_ech = t_fin*length(flute)/max(t);  
flute_morceau = flute_t(t_debut_ech:t_fin_ech); 

% Pisarenko 

[ff_pisa_0, mydsp_pisa_0] = mypisarenko(flute_morceau, 100, fs, 0); 
[ff_pisa_1, mydsp_pisa_1] = mypisarenko(flute_morceau, 100, fs, 1);
[ff_pisa_2, mydsp_pisa_2] = mypisarenko(flute_morceau, 100, fs, 2);
  %%%   - pp : ordre du mod\`ele (choisi de mani\`ere ind\'ependante)
  %%%   - fe : fr\'equence d'\'echantillonnage
  %%%   - rec: 0 : toutes les valeurs propres                           + DSP
  %%%          1 : recursion pour avoir la valeur propre la plus petite + spectre de raies
  %%%          2 : recursion pour avoir la valeur propre la plus petite + enveloppe interpol\'ee

 
% On ne veut prendre que la partie où les fréquences sont positives 
ff_pisa_real_0 = ff_pisa_0(floor(length(ff_pisa_0)/2)+1:end); 
mydsp_pisa_real_0 = mydsp_pisa_0(floor(length(mydsp_pisa_0)/2)+1:end);
figure; 
plot(ff_pisa_real_0, mydsp_pisa_real_0);
title('Paramètre 0') 

% On ne veut prendre que la partie où les fréquences sont positives 
ff_pisa_real_1 = ff_pisa_1(floor(length(ff_pisa_1)/2)+1:end); 
mydsp_pisa_real_1 = mydsp_pisa_1(floor(length(mydsp_pisa_1)/2)+1:end);
figure; 
plot(ff_pisa_real_1, mydsp_pisa_real_1);
title('Paramètre 1') 

% 
ff_pisa_real_2 = ff_pisa_2(floor(length(ff_pisa_2)/2)+1:end); 
mydsp_pisa_real_2 = mydsp_pisa_2(floor(length(mydsp_pisa_2)/2)+1:end);
figure; 
plot(ff_pisa_real_2, mydsp_pisa_real_2);
title('Paramètre 2') 


% On cherche un vecteur de pics de fréquences 
%[pks_pisa_0 loc_pisa_0] = findpeaks(mydsp_pisa_real_0); 
%index_pks_pisa_0 = ff_pisa_real_0(loc_pisa_0)

 
