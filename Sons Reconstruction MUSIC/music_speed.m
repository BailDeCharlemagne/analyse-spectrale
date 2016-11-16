clear all;
close all;

[signal,fe] = audioread('fluteircam.wav');

%%% Setup des constantes du programme
interval_time = .03; % Taille de la fenêtre de calcul de la DSP
time_shift = .01; % Décalage de la fenêtre à chaque itération
df = 0.1; % Précision en fréquence en sortie de mymusic
fmin = 100;
fmax = 800;
tmin = 0.4;
tmax = 1.0;

pic_min_amp = 12 ; % pic_min_amp de l'ordre de 10 donne de bons résultats
signal = signal(fe*tmin + 1:fe*tmax); 

% Génération de la liste des fréquences
LA3 = 440;
freq_list = zeros(1,49);
for k = -24:24
    freq_list(k+25)=LA3 * 2^(k/12);
end


% Préparation de la récupération des pics fréquentiels des DSP mymusic
nPlus = ceil((tmax - tmin)/time_shift);
nMoins = floor(interval_time/time_shift);
numberOfIntervals = nPlus - nMoins +1;
interval_points = floor(interval_time*fe);
n_shift = floor(time_shift*fe);
pics = [];

% elapsed_time = zeros(4,10);
elapsed_time = [];
% Récupération des pics fréquentiels des DSP mymusic

for p = 5:5:50
    for M = 2*p:p:5*p
        tic
        for i = 1 : numberOfIntervals
            signal_iter = signal(1 + (i-1)*n_shift : (i-1)*n_shift + interval_points);
            [ff, mydsp] = mymusic(signal_iter, p, M, fe, df, fmin, fmax);
        end
%         elapsed_time(p/5,M/p-1) = toc;
        elapsed_time = [elapsed_time toc];
    end
end
%%% NOTA %%%
% pics(1,:) = fréquences
% pics(2,:) = amplitudes 
% pics(3,:) = index
% pics(4,:) = temps

% Affichage des temps

figure
plot(elapsed_time,'o');
%     subplot(3,1,1)
%     plot(pics(4,:),pics(1,:),'o')
%     xlabel('Temps (s)');
%     ylabel('Fréquence (Hz)');
%     title(strcat('Pics obtenus avec p=', num2str(p),', M=', num2str(M),' et alpha_lissage = ', num2str(alpha_lissage)));
%     subplot(3,1,2)
%     plot(pics(1,:),log10(pics(2,:)),'o')
%     xlabel('Fréquence (Hz)');
%     ylabel('log(Amplitude)');
%     subplot(3,1,3)
%     plot(pics(4,:),log10(pics(2,:)),'o')
%     xlabel('Temps (s)');
%     ylabel('log(Amplitude)');


