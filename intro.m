function intro

    [xx, fe] = audioread('fluteircam.wav');
    tt = (1:length(xx))/fe;
    plot(tt, xx);
    title('Signal de fl�te');
    xlabel('Time (s)');
    ylabel('Amplitude normalis�e');

end

