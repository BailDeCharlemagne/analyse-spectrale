function intro

    [xx, fe] = audioread('fluteircam.wav');
    tt = (1:length(xx))/fe;
    plot(tt, xx);
    title('Signal de flûte');
    xlabel('Time (s)');
    ylabel('Amplitude normalisée');

end

