[xx,fe] = audioread('fluteircam.wav');

nfreq = 4*64;
decf = 1;
f_Subsampling = 2;

xx2 = xx(1:f_Subsampling:end);
[tfd, t, f] = stft2(xx2, fe, nfreq/f_Subsampling, decf);
f_red = f(nfreq/2-1:end);
abs_TFD = abs(tfd(1:nfreq/f_Subsampling/2+1,:));
imagesc(t,f_red,(flipud(abs_TFD)));