function [toneTxt, tone, duration] = tf_analysis()
    close all;
    clf;

    [xx, fe] = audioread('fluteircam.wav');

    figure(1)
    [xx_reconstructed, fond_F, fond_F_tone, fond_A] = fond_extract(xx, fe, 0.03, 1);
    audiowrite('analyse-spectrale/tf_analysis/tf-fluteircam-reconstruction.wav', xx_reconstructed, fe);
    
    figure(5);
    [tone, duration, toneTxt] = create_partition(fond_F_tone, fond_A, fe*length(fond_F_tone)/length(xx), 1);

    return;
    
    old_xx = xx;
    Pxx = norm(xx)^2/length(xx);
    PxxdB = pow2db(Pxx);
    xx = xx + wgn(length(xx), 1, PxxdB + 15);
    audiowrite('fluteircam-noisy.wav', xx, fe);

    figure(4)
    hold on
    plot(xx)
    plot(old_xx)
    hold off

    figure(2)
    fond_extract(xx, fe, 50, 1);
    figure(3)
    fond_extract(old_xx, fe, 100, 1);
    [xx_reconstructed] = fond_extract(xx, fe, 1000, 1);
    audiowrite('fluteircam-noisy-reconstruction.wav', xx_reconstructed, fe);
end

function [tone, duration, toneTxt] = create_partition(fond_F_tone, fond_A, fe, do_plot)
    fond_F_tone = fond_F_tone(fond_A > 0.01);
    fond_A = fond_A(fond_A > 0.01);
    tone = (fond_F_tone(1));
    duration = (1);
    amplitude = (fond_A(1));
    n = 1;
    for i = 2:length(fond_F_tone)
        if (tone(n) == fond_F_tone(i))
            duration(n) = duration(n) + 1;
            amplitude(n) = ((n-1)*amplitude(n) + fond_A(i)) / n;
        else
            n = n + 1;
            tone(n) = fond_F_tone(i);
            duration(n) = 1;
            amplitude(n) = fond_A(i);
        end
    end
    duration = duration / fe;
    tone = tone(duration > 0.05);
    amplitude = amplitude(duration > 0.05);
    duration = duration(duration > 0.05);
    
    if (do_plot)
        subplot(3,1,1)
        bar(tone, 'blue');
        title('Partition reconstruction');
        ylabel('Difference in semitones');
        subplot(3,1,2)
        bar(duration, 'red');
        ylabel('Duration (s)');
        subplot(3,1,3)
        bar(amplitude, 'green');
        ylabel('Amplitude');
        hold off
    end
    
    tones = cellstr(['La  '; 'Sib '; 'Si  '; 'Do  '; 'Do# '; 'Ré  '; 'Mib '; 'Mi  '; 'Fa  '; 'Fa# '; 'Sol '; 'Sol#']);
    toneTxt = tones(mod(tone, 12)+1);
    %toneTxt = 440*2.^(tone/12);
end

function [xx_reconstructed, fond_F, fond_F_tone, fond_A, abs_TFD] = fond_extract(xx, fe, tw, do_plot)

    f_Subsampling = 4;   % sub-sampling factor
    nfreq = 64*64;       % number of samples to compute in frequency
    decf = 16;            % sub-sampling factor in time of the stft

    % Computes the short-time Fourier transform
    xx2 = xx(1:f_Subsampling:end);
    w = round(tw * (fe/f_Subsampling))+1;
    [tfd, t, f, w2] = stft2(xx2, fe/f_Subsampling, nfreq, decf, ones(w,1));
    f_red = f(nfreq/2:end);
    abs_TFD = abs(tfd(1:nfreq/2+1,:));
    
    % Extract the fundamental tone and amplitude
    [fond_A, fond_F] = max(flipud(abs_TFD));
    fond_F = f_red(fond_F);
    fond_F_tone = round(log2(fond_F/440)*12);
%     figure(6)
%     hold on
%     for t = -30:30
%         toto = fond_F(fond_F_tone==t);
%         histogram(toto);
%     end
%     hold off
% 
%     figure(1)
    fond_F = 440*2.^(fond_F_tone/12);
    
    do_plot = 1
    if (do_plot)
%         subplot(1,1,1);
%         imagesc(t, f_red, flipud(abs_TFD));
%         str = sprintf('Short-time Fourrier transform (w = %ims)', round(w/(fe/f_Subsampling)*1000));
%         title(str);
%         xlabel('Time (s)');
%         ylabel('Frequency (Hz)');

        subplot(2,1,1);
        title('Fundamental tone extraction');
        plot(t, fond_F_tone);
        ylabel('Difference in semitones with LA 440');
        xlabel('Time (s)');
        
        subplot(2,1,2);
        plot(t, fond_A);
        ylabel('Amplitude');
        xlabel('Time (s)');
    end

    % Reconstructs the signal
    xx_reconstructed = zeros(length(xx), 1);
    display(size(fond_F));
    for ii = 1:length(xx)
        %tt = ii / fe;
        %ii2 = round(tt * f_Subsampling * decf * 8);
        ii2 = round(ii * length(fond_A) / length(xx));
        ii2 = max(min(ii2, length(fond_A)), 1);
        xx_reconstructed(ii,1) = fond_A(ii2) * sin (2*pi*fond_F(ii2) * ii/fe);
    end
    
    figure()
    hold on
    plot((1:length(xx))/fe, xx_reconstructed)
    plot((1:length(xx))/fe, xx)
    legend('Signal reconstruit', 'Signal initial')
    hold off
    
    % Correct the signal power
    C = (norm(xx) / norm(xx_reconstructed))^2;
    xx_reconstructed = xx_reconstructed * C;
    fond_A = fond_A * C;
end
