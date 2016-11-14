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
%figure;
%plot(t, flute);
%title('Signal pur de la flute - non bruité'); 

 
####################################################################
##########       PISARENKO SUR LE SIGNAL ENTIER       ##############
####################################################################

close all; 
order = [1 2 10 50 100];
parameter = 1 ; % 0 donne la DSP, 1 le spectre de raies, 2 l'enveloppe interpolée
plot_or_not_plot = 0;

% Color Map for plots
G = fliplr(linspace(0,1,length(order))) .';
myGmap = horzcat(zeros(size(G)), G, zeros(size(G))); 
 
figure;
power_density = []; 

for i = 1:length(order)
  % Utilisation de Pisarenko 
  [ff_pisa, densite_puissance] = mypisarenko(flute_t, order(i), fs, 0); % donne ka densité spec
  [ff_pisa, spectre_pisa] = mypisarenko(flute_t, order(i), fs, 1); % donne spectre de raies
  [ff_pisa, enveloppe_interpolee] = mypisarenko(flute_t, order(i), fs, 2); % donne enveloppe interpolée
  subplot(length(order),1,i);
  plot(ff_pisa, spectre_pisa, 'color', myGmap(i,:)); axis tight; hold on;
  plot(ff_pisa, enveloppe_interpolee, 'color', myGmap(i,:));
  hold off;
  str = ['Spectre de raies et Enveloppe Interpolée pour p = ', num2str(order(i))];
  % Quelle est la densité spectrale de puissance ?
  power_density = [power_density densite_puissance]; 
  % Quelles sont les fréquences pour lesquelles on observe des pics ?
  [pks_pisa loc_pks] = findpeaks(spectre_pisa);
  %peaks{i} = (pks_pisa,ff_pisa(loc_pks));
  M = zeros(length(pks_pisa), 3); % matrix that will contain the location of the peaks, their frequency, and the amplitude
  M(:,1) = pks_pisa; M(:,2) = loc_pks; M(:,3) = ff_pisa(loc_pks);
  N = sortrows(M, -1); 
  freq{i} = N(:,3)(N(:,3)>0); 
  title(str);
  xlim([-5000, 5000]); 
end;

freq
power_density