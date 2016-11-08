[xx,fe] = audioread('fluteircam.wav');

f_Subsampling = 4;   % sub-sampling factor
nfreq = 64*64;       % number of samples to compute in frequency
decf = 32;           % sub-sampling factor in time of the stft
w = 500;             % length of the window

% Computes the short-time Fourier transform
xx2 = xx(1:f_Subsampling:end);
[tfd, t, f] = stft2(xx2, fe/f_Subsampling, nfreq, decf, w);
f_red = f(nfreq/2+2:end);
abs_TFD = abs(tfd(1:nfreq/2+1,:));

% Extract the fundamental tone and amplitude
[fond_A, fond_F] = max(flipud(abs_TFD));
fond_F = f_red(fond_F);
fond_F_tone = round(log2(fond_F/440)*12);
fond_F = 440*2.^(fond_F_tone/12);

% Plots
close all;
clf;

figure(1)
subplot(2,1,1);
imagesc(t,f_red,(flipud(abs_TFD)));
title('Short-time Fourrier transform');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(2,1,2);
hold on
title('Fundamental tone extraction');
yyaxis left
plot(t, fond_F_tone);
ylabel('Difference in semitones with LA 440');
yyaxis right
plot(t, fond_A);
ylabel('Amplitude');
hold off
xlabel('Time (s)');

% Reconstructs the signal with the fondamental only
figure(2)
hold on
plot((1:length(xx))*fe, xx)
for ii = 1:length(xx)
    tt = ii / fe;
    tt_red = round(tt * f_Subsampling * decf * 2);
    tt_red = max(1,min(tt_red, length(fond_A)));
    xx(ii,1) = fond_A(tt_red) * sin (2*pi*fond_F(tt_red) * tt);
end
plot((1:length(xx))*fe, xx)
hold off
xlabel('Time (s)');
ylabel('Amplitude');
legend('Raw signal', 'Reconstructed signal');
audiowrite('fluteircam-reconstruction.wav', xx, fe);