clear all;
close all;

[signal,fe] = audioread('fluteircam.wav');

%%% Setup des options du programme (souvent modifiées)
p = 220;
M = 300; % Prendre M >> p
alpha_lissage = 0;
interpolation = 1;
write_fondamental = 0;
write_rounded_fondamental = 0;

display_all_peaks = 1;
display_fondamental = 1;
display_rounded_fondamental = 1;

%%% Setup des constantes du programme
interval_time = .03; % Taille de la fenêtre de calcul de la DSP
time_shift = .01; % Décalage de la fenêtre à chaque itération
df = 0.1; % Précision en fréquence en sortie de mymusic
fmin = 100;
fmax = 800;
tmin = 4;
tmax = 6;

pi = 3.1415;
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

% Récupération des pics fréquentiels des DSP mymusic
for i = 1 : numberOfIntervals
    signal_iter = signal(1 + (i-1)*n_shift : (i-1)*n_shift + interval_points);
    [ff, mydsp] = mymusic(signal_iter, p, M, fe, df, fmin, fmax);
    [val,locs] = findpeaks(mydsp);
    f_pks = [fmin + (locs-1)*df];
    vect_i = ones(length(val),1)*i;
    vect_t = ones(length(val),1)*(tmin + interval_time * 0.5 + time_shift * (i-1));
    pics = [pics [f_pks, val,vect_i,vect_t]'];
end 
%%% NOTA %%%
% pics(1,:) = fréquences
% pics(2,:) = amplitudes 
% pics(3,:) = index
% pics(4,:) = temps

% Suppression des pics de trop faible amplitude
suppr_pics = [];
for i = 1:length(pics)
    if pics(2,i) < pic_min_amp
       suppr_pics = [suppr_pics i];
    end
end
pics(:,suppr_pics) = [];

% Obtention de l'évolution temporelle de la fréquence fondamentale
fondamental = [];
iterator = 1;
index = 1;
while( index <= numberOfIntervals && iterator <= size(pics,2))
    iterator_fondamental = iterator;
    while(iterator <= size(pics,2) && pics(3,iterator) == index);
        if pics(2,iterator) > pics(2,iterator_fondamental)
            iterator_fondamental = iterator;
        end
        iterator = iterator + 1 ;
    end
    fondamental = [fondamental pics(:,iterator_fondamental)];
    index = index+1;
end

% vect_440 = ones(size(fondamental(1,:))) * 440
% vect_mean_fond = ones(size(fondamental(1,:))) * mean(fondamental(1,:))
% figure
% plot(fondamental(4,:),fondamental(1,:), 'o')
% hold on
% plot(fondamental(4,:),vect_440, 'r')
% hold on
% plot(fondamental(4,:),vect_mean_fond, 'g')
% xlabel('Temps (s)');
% ylabel('Fréquence (Hz)');
%  title(strcat('Fréquences fondamentales obtenues avec p=', num2str(p),', M=', num2str(M)));

% Projection de la fondamentale sur le demi-tons
rounded_fondamental = fondamental(:,:);
for iter = 1:size(rounded_fondamental,2)
    [m, index_list] = min(abs(freq_list-rounded_fondamental(1,iter)));
    rounded_fondamental(1,iter) = freq_list(index_list);
    rounded_fondamental(2,iter) = 1000;
end

% Obtention du biais et de la variance dans les données
bias = mean(abs(440-fondamental(1,:)));
variance = var(fondamental(1,:));

% Lissage en amplitude et en fréquence via le coefficient alpha_lissage
if (alpha_lissage>0)
    for i= 2 : size(fondamental,2)
        fondamental(2,i) = alpha_lissage * fondamental(2,i-1) + (1-alpha_lissage) * fondamental(2,i);
        fondamental(1,i) = alpha_lissage * fondamental(1,i-1) + (1-alpha_lissage) * fondamental(1,i);
    end
end

% Affichage de tous les pics
if(display_all_peaks==1)
    figure
    subplot(3,1,1)
    plot(pics(4,:),pics(1,:),'o')
    xlabel('Temps (s)');
    ylabel('Fréquence (Hz)');
    title(strcat('Pics obtenus avec p=', num2str(p),', M=', num2str(M),' et alpha_lissage = ', num2str(alpha_lissage)));
    subplot(3,1,2)
    plot(pics(1,:),20*log10(pics(2,:)),'o')
    xlabel('Fréquence (Hz)');
    ylabel('Puissance (dB)');
    subplot(3,1,3)
    plot(pics(4,:),20*log10(pics(2,:)),'o')
    xlabel('Temps (s)');
    ylabel('Puissance (dB)');
end


% Affichage du fondamental
if(display_fondamental==1)
    figure
    subplot(3,1,1)
    plot(fondamental(4,:),fondamental(1,:))
    xlabel('Temps (s)');
    ylabel('Fréquence (Hz)');
    title(strcat('Fondamental obtenu avec p=', num2str(p),', M=', num2str(M),' et alpha_lissage = ', num2str(alpha_lissage)));
    subplot(3,1,2)
    plot(fondamental(1,:),20*log10(fondamental(2,:)),'o')
    xlabel('Fréquence (Hz)');
    ylabel('Puissance lissée(dB)');
    subplot(3,1,3)
    plot(fondamental(4,:),20*log10(fondamental(2,:)))
    xlabel('Temps (s)');
    ylabel('Puissance lissée(dB)');
end

% Affichage du fondamental arrondi
if(display_rounded_fondamental==1)
    figure
    subplot(3,1,1)
    plot(rounded_fondamental(4,:),rounded_fondamental(1,:))
    xlabel('Temps (s)');
    ylabel('Fréquence (Hz)');
    title(strcat('Fondamental arrondi au demi-ton le plus proche, Résultats obtenus avec p=', num2str(p),' et M=', num2str(M),' et alpha_lissage = ', num2str(alpha_lissage)));
    subplot(3,1,2)
    plot(rounded_fondamental(1,:),20*log10(rounded_fondamental(2,:)),'o')
    xlabel('Fréquence (Hz)');
    ylabel('Puissance lissée (dB)');
    subplot(3,1,3)
    plot(rounded_fondamental(4,:),20*log10(rounded_fondamental(2,:)))
    xlabel('Temps (s)');
    ylabel('Puissance lissée (dB)');
end

% Ecriture du fondamental sans interpolation dans un fichier .wav
if (interpolation == 0 && write_fondamental == 1)
    xx = zeros(1, numberOfIntervals * n_shift);
    phi=0;
    for it = 1 : numberOfIntervals
        phi = phi + (2*pi*fondamental(1,it)*n_shift)/fe;
        for i = 1:n_shift    
            xx(n_shift * (it-1) + i) = sqrt(fondamental(2,it))* cos(2*pi*fondamental(1,it)*i/fe+phi);
        end
    end
    % strcat('fluteircam-MUSIC_reconstruction_tmin=',num2str(tmin),'_tmax=',num2str(tmax),'_p=',num2str(p),'_M=',num2str(M),'_alpha=',num2str(alpha_lissage),'.wav');
    audiowrite('fluteircam-MUSIC_reconstruction_tmin=_tmax=_p=_M=_alpha_lissage=.wav', xx, fe);
end

% % % % Ecriture du fondamental avec interpolation dans un fichier .wav
% % % if(interpolation == 1 && write_fondamental == 1)
% % %     xx = zeros(1, numberOfIntervals * n_shift);
% % %     phi=0;
% % %     for it = 1 : numberOfIntervals-1
% % %         phi = phi + (2*pi*fondamental(1,it)*n_shift)/fe;
% % %         for i = 1:n_shift    
% % %             amplitude = (n_shift-i)/n_shift * sqrt(fondamental(2,it)) + i/n_shift * sqrt(fondamental(2,it+1));
% % %             xx(n_shift * (it-1) + i) = amplitude * cos(2*pi*fondamental(1,it)*i/fe+phi);
% % %         end
% % %     end
% % %     it = numberOfIntervals;
% % %     phi = phi + (2*pi*fondamental(1,it)*n_shift)/fe;
% % %     for i = 1:n_shift    
% % %         amplitude = sqrt(fondamental(2,it));
% % %         xx(n_shift * (it-1) + i) = amplitude * cos(2*pi*fondamental(1,it)*i/fe+phi);
% % %     end
% % %     % strcat('fluteircam-MUSIC_reconstruction_interpolation_tmin=',num2str(tmin),'_tmax=',num2str(tmax),'_p=',num2str(p),'_M=',num2str(M),'_alpha=',num2str(alpha_lissage),'.wav');
% % %     audiowrite('fluteircam-MUSIC_reconstruction_interpolation_tmin=_tmax=_p=_M=_alpha_lissage=.wav', xx, fe);
% % % end

% Ecriture du fondamental arrondi sans interpolation dans un fichier .wav
if (interpolation == 0 && write_rounded_fondamental == 1)
    xx = zeros(1, numberOfIntervals * n_shift);
    phi=0;
    for it = 1 : numberOfIntervals
        phi = phi + (2*pi*rounded_fondamental(1,it)*n_shift)/fe;
        for i = 1:n_shift    
            xx(n_shift * (it-1) + i) = sqrt(rounded_fondamental(2,it))* cos(2*pi*rounded_fondamental(1,it)*i/fe+phi);
        end
    end
    % strcat('tmin=',num2str(floor(tmin)),'_tmax=',num2str(ceil(tmax)),'_p=',num2str(p),'_M=',num2str(M),'_alpha=',num2str(alpha_lissage),'.wav');
    audiowrite('fluteircam-MUSIC_rounded_reconstruction_tmin=_tmax=_p=_M=_alpha_lissage=.wav', xx, fe);
end

% Ecriture du fondamental arrondi avec interpolation dans un fichier .wav
if(interpolation == 1 && write_rounded_fondamental == 1)
    xx = zeros(1, numberOfIntervals * n_shift);
    phi=0;
    for it = 1 : numberOfIntervals-1
        phi = phi + (2*pi*rounded_fondamental(1,it)*n_shift)/fe;
        for i = 1:n_shift    
            amplitude = (n_shift-i)/n_shift * sqrt(rounded_fondamental(2,it)) + i/n_shift * sqrt(rounded_fondamental(2,it+1));
            xx(n_shift * (it-1) + i) = amplitude * cos(2*pi*rounded_fondamental(1,it)*i/fe+phi);
        end
    end
    it = numberOfIntervals;
    phi = phi + (2*pi*rounded_fondamental(1,it)*n_shift)/fe;
    for i = 1:n_shift    
        amplitude = sqrt(rounded_fondamental(2,it));
        xx(n_shift * (it-1) + i) = amplitude * cos(2*pi*rounded_fondamental(1,it)*i/fe+phi);
    end
    % strcat('fluteircam-MUSIC_reconstruction_interpolation_tmin=',num2str(tmin),'_tmax=',num2str(tmax),'_p=',num2str(p),'_M=',num2str(M),'_alpha=',num2str(alpha_lissage),'.wav');
%     audiowrite('fluteircam-MUSIC_rounded_reconstruction_interpolation_tmin=_tmax=_p=_M=_alpha_lissage=.wav', xx, fe);
temps_xx = (1:size(xx,2))/fe+0.1;
temps_signal = (1:size(signal,2))/fe;

figure
subplot(2,1,1)
plot(temps_xx,xx, 'g')
hold on
plot(temps_signal,signal, 'r')
xlabel('Temps (s)');
ylabel('Signal');
legend('signal reconstruit','signal original')
title(strcat('Reconstruction obtenue avec p=', num2str(p),', M=', num2str(M), 'en projetant sur les demi-tons les plus proches, à puissance constante' ));
end




% Ecriture du fondamental avec interpolation dans un fichier .wav
if(interpolation == 1 && write_fondamental == 1)
    xx = zeros(1, numberOfIntervals * n_shift);
    phi=0;
    for it = 1 : numberOfIntervals-1
        phi = phi + (2*pi*fondamental(1,it)*n_shift)/fe;
        for i = 1:n_shift    
            amplitude = ((n_shift-i)/n_shift * sqrt(fondamental(2,it)) + i/n_shift * sqrt(fondamental(2,it+1)))/500;
%             amplitude = ((n_shift-i)/n_shift * sqrt(fondamental(2,it)) + i/n_shift * sqrt(fondamental(2,it+1)));
            frequence = (n_shift-i)/n_shift * fondamental(1,it) + i/n_shift * fondamental(1,it+1);
            xx(n_shift * (it-1) + i) = amplitude * cos(2*pi*frequence*i/fe+phi);
        end
    end
    it = numberOfIntervals;
    phi = phi + (2*pi*fondamental(1,it)*n_shift)/fe;
    for i = 1:n_shift    
        amplitude = sqrt(fondamental(2,it))/500;
%         amplitude = sqrt(fondamental(2,it));
        frequence =  fondamental(1,it);
        xx(n_shift * (it-1) + i) = amplitude * cos(2*pi*fondamental(1,it)*i/fe+phi);
    end
    % strcat('fluteircam-MUSIC_reconstruction_interpolation_tmin=',num2str(tmin),'_tmax=',num2str(tmax),'_p=',num2str(p),'_M=',num2str(M),'_alpha=',num2str(alpha_lissage),'.wav');
%     audiowrite('fluteircam-MUSIC_reconstruction_interpolation_frequence_tmin=_tmax=_p=_M=_alpha_lissage=.wav', xx, fe);
% temps_xx = (1:size(xx,2))/fe+0.1;
% temps_signal = (1:size(signal,2))/fe;
% 
% subplot(2,1,2)
% plot(temps_xx,xx, 'g')
% hold on
% plot(temps_signal,signal, 'r')
% xlabel('Temps (s)');
% ylabel('Signal');
% legend('signal reconstruit','signal original')
% title(strcat('Reconstruction obtenue avec p=', num2str(p),', M=', num2str(M), ' avec interpolation, lissage et diminution d amplitude' ));
end



