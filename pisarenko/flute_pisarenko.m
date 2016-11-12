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
%figure; 
%plot(t, noisy_flute) 
%title('Signal bruité de la flute');  


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

%%% Pisarenko sur le signal entre 4 et 6 secondes

# Pisarenko 0
[ff_pisa0, dsp] = pisarenko_real(flute_morceau, 100, fs, 0, 0); % donne densité de puissance spectrale

# Pisarenko 1 
[ff_pisa1, mydsp_pisa1] = pisarenko_real(flute_morceau, 100, fs, 1, 0); % donne spectre de raies
% On cherche un vecteur de pics de fréquences 
[pks_pisa1 loc_pisa1] = findpeaks(mydsp_pisa1); 
index_pks_pisa1 = ff_pisa1(loc_pisa1); 

# Pisarenko 2 : 
[ff_pisa2, mydsp_pisa2] = pisarenko_real(flute_morceau, 100, fs, 2,0); % donne enveloppe interpolée
% On cherche un vecteur de pics de fréquences 
[pks_pisa2 loc_pisa2] = findpeaks(mydsp_pisa2); 
index_pks_pisa_2 = ff_pisa2(loc_pisa2);

% Cette figure va nous donner le spectre de raies et l'enveloppe interpolée
figure;
plot(ff_pisa1, mydsp_pisa1, mydsp_pisa2);  

 
###########################################################################
##########       PISARENKO SUR LES PETITS INTERVALLES        ##############
###########################################################################

close all; 

%% Construction des petits signaux 

interval = 30*10^(-3); 
interval_ech = interval*fs;
nbIntervals = floor(length(flute_t)/interval_ech); 
flute_ech = zeros(nbIntervals, interval_ech);
for i=1:nbIntervals
  flute_ech(i,:) = flute_t((i-1)*interval_ech+1:i*interval_ech); 
end; 

%%% Pisarenko sur ces petits signaux 

%% Matrice pour construire le temps fréquence amplitude
M = zeros(nbIntervals, fs/8);

for i=1:nbIntervals
  # Pisarenko 1 
  [ff_pisa, mydsp_pisa] = pisarenko_real(flute_ech(i,:), 100, fs, 1, 0); % donne spectre de raies
  mydsp_pisa1 = zeros(1,mydsp_pisa);
  for k=1:length(mydsp_pisa) 
    if (mydsp_pisa(k) > 10^(-7))
      mydsp_pisa1(k) = mydsp_pisa(k);
    end; 
  end; 
  M(i,:) = transpose(mydsp_pisa1(1:fs/8)); 
  str = ['done ', num2str(i)];
  disp(str)
end;

figure;
mesh(M); 
xlabel('Fréquence'); ylabel('Intervalle n°'); zlabel('Amplitude de la DSP'); 
rotate3d on; 


%% On cherche un vecteur de pics de fréquences 
%[pks_pisa1 loc_pisa1] = findpeaks(mydsp_pisa1); 
%index_pks_pisa1 = ff_pisa1(loc_pisa1); 
%
%# Pisarenko 2 : 
%[ff_pisa2(i,:), mydsp_pisa2(i,:)] = pisarenko_real(flute_ech(i,:), 100, fs, 2,0); % donne enveloppe interpolée
%% On cherche un vecteur de pics de fréquences 
%[pks_pisa2 loc_pisa2] = findpeaks(mydsp_pisa2); 
%index_pks_pisa_2 = ff_pisa2(loc_pisa2);
%
%figure;
%plot(ff_pisa1, mydsp_pisa1, mydsp_pisa2);  
