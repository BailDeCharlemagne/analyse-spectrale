%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PREPROCESSING %%%%%%%%
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


###########################################################################
##########            PISARENKO SUR 0.6 A 1.4 SECONDES          ###########
###########################################################################


close all; 

% données temporelles
t_debut = 0.4; %en sec 
t_fin = 1; %en sec
saut = 10*10^(-3); %en sec
intervalle = 30*10^(-3); %en sec

% données échantillonnées 
t_debut_ech = t_debut*length(flute)/max(t); 
t_fin_ech = t_fin*length(flute)/max(t);
saut_ech = saut*length(flute)/max(t); 
intervalle_ech = intervalle*length(flute)/max(t);
nb_ech = (t_fin - t_debut)/saut;

% morceau de flute auquel on s'intéresse 
flute_morceau = flute_t(t_debut_ech:t_fin_ech); 
flute_ech = ones(nb_ech, intervalle_ech);

% On splite sur chaque échantillon 
for i = 1:int8(nb_ech)
  flute_ech(i,:) = flute_t(1,t_debut_ech+saut_ech*(i-1):t_debut_ech+saut_ech*(i-1)+(intervalle_ech-1));
end;  



%%% Pisarenko sur les signaux

order = 20; 
parameter = 1 ; % 0 donne la DSP, 1 le spectre de raies, 2 l'enveloppe interpolée
plot_or_not_plot = 0;

% Matrice pour construire le temps fréquence amplitude
M_full = zeros(int8(nb_ech), fs/8); %lignes n° échantillon + colonnes fréquence 

tic();

for i=1:int8(nb_ech)
  [ff_pisa, mydsp_pisa] = pisarenko_real(flute_ech(i,:), order, fs, parameter, plot_or_not_plot); % donne spectre de raies
  M_full(i,:) = transpose(mydsp_pisa(1:fs/8)); 
end;

[max_pics, max_indices] = sort(M_full, 2, 'descend');
max_pics = max_pics(:, 1);
max_indices = max_indices(:, 1);
freq_max_ampl = ff_pisa(max_indices(:,1)); 

elapsed_time = toc();

biais = 0; 
variance = 0; 
freq_th = 440;

for i=1:nb_ech
  biais += abs(freq_max_ampl(i)-freq_th);
  variance += (freq_max_ampl(i)-freq_th)^2; 
end; 

biais = biais/nb_ech;
variance = sqrt(1/nb_ech * variance);

order
elapsed_time
moyenne = mean(freq_max_ampl) 
biais 
variance
