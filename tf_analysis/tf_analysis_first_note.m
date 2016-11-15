function [bias, variance] = tf_analysis_first_note()
    close all;
    clf;

    [xx, fe] = audioread('fluteircam.wav');

    ti = 0.4;             % start (s)
    tf = 1.0;             % end (s)
    f_Subsampling = 1;    % sub-sampling factor
    nfreq = 2048*64;       % number of samples to compute in frequency
    decf = 32;            % sub-sampling factor in time of the stft
    tw = 0.03;

    % Computes the short-time Fourier transform
    xx2 = xx(ti*fe:f_Subsampling:tf*fe);
    w = round(tw * (fe/f_Subsampling))+1;
    [tfd, t, f] = stft2(xx2, fe/f_Subsampling, nfreq, decf, w);
    f_red = f(nfreq/2:end);
    abs_TFD = abs(tfd(1:nfreq/2+1,:));

    [fond_A, fond_F] = max(flipud(abs_TFD));
    fond_F = f_red(fond_F);
    
    figure(1)
    
    subplot(1,2,1)
    imagesc(t+ti, f_red, flipud(abs_TFD));
    str = sprintf('Short-time Fourrier transform (w = %ims)', round(w/(fe/f_Subsampling)*1000));
    title(str);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    subplot(1,2,2)
    hold on
    plot(t+ti, fond_F)
    plot(t+ti, ones(size(t))*440)
    hold off
    xlabel('Time (s)');
    ylabel('Fundamental frequency (Hz)');
    
    bias = abs(440-mean(fond_F))
    variance = var(fond_F)
end
