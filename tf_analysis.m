function tf_analysis()
    close all;
    clf;

    [xx, fe] = audioread('fluteircam.wav');

    figure(1)
    [xx_reconstructed] = fond_extract(xx, fe, 500, 1);
    audiowrite('fluteircam-reconstruction.wav', xx_reconstructed, fe);

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

function [xx_reconstructed, fond_F, fond_F_tone, fond_A, abs_TFD] = fond_extract(xx, fe, w, do_plot)

    f_Subsampling = 4;   % sub-sampling factor
    nfreq = 64*64;       % number of samples to compute in frequency
    decf = 16;           % sub-sampling factor in time of the stft

    % Computes the short-time Fourier transform
    xx2 = xx(1:f_Subsampling:end);
    [tfd, t, f] = stft2(xx2, fe/f_Subsampling, nfreq, decf, w);
    f_red = f(nfreq/2:end);
    abs_TFD = abs(tfd(1:nfreq/2+1,:));

    % Extract the fundamental tone and amplitude
    [fond_A, fond_F] = max(flipud(abs_TFD));
    fond_F = f_red(fond_F);
    fond_F_tone = round(log2(fond_F/440)*12);
    fond_F = 440*2.^(fond_F_tone/12);
    
    if (do_plot)
        subplot(2,1,1);
        imagesc(t, f_red, flipud(abs_TFD));
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
    end

    % Reconstructs the signal
    xx_reconstructed = zeros(length(xx), 1);
    for ii = 1:length(xx)
        tt = ii / fe;
        tt_red = round(tt * f_Subsampling * decf * 8);
        tt_red = max(1,min(tt_red, length(fond_A)));
        xx_reconstructed(ii,1) = fond_A(tt_red) * sin (2*pi*fond_F(tt_red) * tt);
    end
    
    % Correct the signal power
    C = (norm(xx) / norm(xx_reconstructed))^2;
    xx_reconstructed = xx_reconstructed * C;
    fond_A = fond_A * C;
end
