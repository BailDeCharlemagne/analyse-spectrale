%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PREPROCESSING %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
 
%% Partition de référence 
partition_ref = [440,494,523,659,415,440,349,330,370,415,440,494,587,698,659,587,494,523,440,494] 

%% Gamme 
n = -10:1:20;
f0 = 440; % La 
gamme = [];
for i=1:length(n)
  gamme(i) = f0*2^(n(i)/12)
end; 
 
 
###########################################################################
##########       PISARENKO SUR LES PETITS INTERVALLES        ##############
###########################################################################

close all; 
order = 70; 
parameter = 1 ; % 0 donne la DSP, 1 le spectre de raies, 2 l'enveloppe interpolée
plot_or_not_plot = 0;


%% Construction des petits signaux 

interval = 30*10^(-3); 
interval_ech = interval*fs;
nbIntervals = floor(length(flute_t)/interval_ech); 
flute_ech = zeros(nbIntervals, interval_ech);
for i=1:nbIntervals
  flute_ech(i,:) = flute_t((i-1)*interval_ech+1:i*interval_ech); 
end; 

%% Pisarenko sur ces petits signaux 

% Matrice pour construire le temps fréquence amplitude
M_full = zeros(nbIntervals, fs/8);

for i=1:nbIntervals
  # Pisarenko 1 
  [ff_pisa, mydsp_pisa] = pisarenko_real(flute_ech(i,:), order, fs, parameter, plot_or_not_plot); % donne spectre de raies
  mydsp_pisa1 = zeros(1,mydsp_pisa);
  for k=1:length(mydsp_pisa) 
    if (mydsp_pisa(k) > 10^(-7))
      mydsp_pisa1(k) = mydsp_pisa(k);
    else mydsp_pisa1(k) = 0;
    end; 
  end; 
  M_full(i,:) = transpose(mydsp_pisa1(1:fs/8)); 
  str = ['done ', num2str(i)];
  disp(str)
end;

% Récupérer la fréquence dont l'amplitude est max pour chaque intervalle

close all;

[max_pics, max_indices] = sort(M_full, 2, 'descend');
max_pics = max_pics(:, 1:4);
max_indices = max_indices(:, 1:4);
freq_max_ampl = ff_pisa(max_indices); %pour chaque nbInterval, 4 fréquences d'amplitude max
freq_fund= freq_max_ampl(:,1);

%% Traitement des résultats pour display le graphe Temps-fréquence 
figure; 
title('Temps-fréquence du signal - ordre 50');

% Signal brut 
subplot(3,1,1)
plot((1:nbIntervals)*interval,freq_fund, 'b');
ylim([0 1000]); 
legend('Signal non-moyenne non-filtre');
xlabel('Temps (s)'); ylabel('Fréquence (Hz)');

% Arrondi au demi-ton le plus proche 
n(1:nbIntervals) = round(12*log2(freq_fund/f0)); 
freq_fund_round(1:nbIntervals) = f0*2.^(n(1:nbIntervals)/12);
subplot(3,1,2);
plot((1:nbIntervals)*interval,freq_fund_round, 'r');
ylim([0 1000]); 
legend('Signal approché au demi-ton le plus proche'); 
xlabel('Temps (s)'); ylabel('Fréquence (Hz)');

% Filtre pour les fréquences trop hautes  
freq_max = 800; 
freq_fund_round_filter = zeros(size(freq_fund_round));
freq_fund_round_filter(1) = freq_fund_round(1);  
for i=2:length(freq_fund_round)
  if freq_fund_round(i) > freq_max 
      freq_fund_round_filter(i) = freq_fund_round_filter(i-1);
  else
    freq_fund_round_filter(i) = freq_fund_round(i); 
  end;
end; 
subplot(3,1,3);
plot((1:nbIntervals)*interval,freq_fund_round_filter);
ylim([0 1000]); 
legend('Signal approché au demi-ton, filtré sur les HF'); 
xlabel('Temps (s)'); ylabel('Fréquence (Hz)');
title('Temps-fréquence par Pisarenko - Ordre 50'); 

%% Display pour slide 17
%
%figure;
%mesh(M_full); 
%xlabel('Fréquence'); ylabel('Intervalle n°'); zlabel('Amplitude de la DSP'); 
%rotate3d on; 

%% Display en imagesc pour slide 17
%
%figure; 
%imagesc(transpose(M_full));
%xlabel('Interval n°'); ylabel('Fréquence (Hz)')
%title('Temps-Fréquence pour Pisarenko - ordre 50')


###########################################################################
##########                    RECONSTRUCTION                 ##############
###########################################################################