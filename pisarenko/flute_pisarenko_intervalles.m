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
  

 
###########################################################################
##########       PISARENKO SUR LES PETITS INTERVALLES        ##############
###########################################################################

close all; 
order = 1; 
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

% Pisarenko sur ces petits signaux 

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

figure; 
imagesc(transpose(M_full));

% TO DO 
A = [1 2 3 4 8 5 2 9 ; 5 6 4 7 9 1 1 2]
[maxvals, maxindices] = sort(A, 2, 'descend');
maxvals = maxvals(:, 1:4)
maxindices = maxindices(:, 1:4)
%


%% Display pour obtenir slide 17
%
%figure;
%mesh(M_full); 
%xlabel('Fréquence'); ylabel('Intervalle n°'); zlabel('Amplitude de la DSP'); 
%rotate3d on; 

%% Matrice pour construire le temps fréquence amplitude
%M_relevant = zeros(nbIntervals, fs/8);
%
%for i=1:nbIntervals
%  # Pisarenko 1 
%  [ff_pisa, mydsp_pisa] = pisarenko_real(flute_ech(i,:), order, fs, parameter, plot_or_not_plot); % donne spectre de raies
%  mydsp_pisa(mydsp_pisa < 10^-9) = NaN;
%  M_relevant(i,:) = transpose(mydsp_pisa(1:fs/8)); 
%  str = ['done ', num2str(i)];
%  disp(str)
%end;
%
%% Display
%
%figure;
%mesh(M_relevant); 
%xlabel('Fréquence'); ylabel('Intervalle n°'); zlabel('Amplitude de la DSP'); 
%rotate3d on; 

