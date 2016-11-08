% Add Noise to Flute

Eb_No = 0.05; %en dB 
Eb = (norm(flute)^2)/length(flute);
No = Eb_No*log10(Eb); 
variance = No/2; %Revoir avec Telecom
bruit = randn(length(flute),1)*variance; 
noisy_flute = bruit + flute;
figure; 
plot(t, noisy_flute) 
title('Signal bruité de la flute');  


%% FFT du signal 

figure;
hold on; 

noisy_fft_complex = abs(fft(noisy_flute)); 
noisy_fft_real = noisy_fft_complex(1:length(noisy_flute)/2); 
freq_red = fs*(0:length(noisy_flute)/2-1)/length(noisy_flute);
plot(freq_red, noisy_fft_real, 'r');

y_fft_complex = abs(fft(flute_t)); 
y_fft_real = y_fft_complex(1:length(flute)/2); 
freq_red = fs*(0:length(flute)/2-1)/length(flute);
plot(freq_red, y_fft_real, 'b');
hold off;

audiowrite('fluteircam-bruité.wav', noisy_flute, fs);
