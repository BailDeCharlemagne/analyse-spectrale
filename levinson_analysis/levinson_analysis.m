function [fmax, amax, pmax, tmax] = levinson_analysis()

    clear all;
    close all;
    clf;

    [xx, fe] = audioread('fluteircam.wav');
    
    ti = 3.8; % seconds
    tf = 6.0; % seconds
    pp = 150;
    
    %xx2 = flipud(xx(round(ti*fe):round(tf*fe)));
    xx2 = flipud(xx);
    
    figure(1)
    hold on
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx2, 1, fe, 1);
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx2, 5, fe, 1);
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx2, 10, fe, 1);
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx2, 50, fe, 1);
    [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx2, 150, fe, 1);
    legend('1', '5', '10', '50', '150');
    hold off
    
    dt = 0.03;
    L = floor(dt*fe);
    fmax = [];
    amax = [];
    pmax = [];
    tmax = [];
    for ii = 1:L:(length(xx2)-3*L)
        xx3 = xx2(ii:ii+3*L);
        [aa, sigma2, ref, ff, mydsp] = mylevinsondurbin(xx3, pp, fe, 0);
        [peaks locs] = findpeaks(mydsp);
        
        jj = floor(ii/L)+1;
        fmax(jj) = locs(1) * (ff(2)-ff(1));
        amax(jj) = peaks(1);
        pmax(jj) = norm(xx3)^2 / length(xx3);
        tmax(jj) = ii/fe;

        wb = waitbar(ii/length(xx2));
    end
    delete(wb);

    figure(2);
    title('Vibrato analysis (Levinson150)');
    
    subplot(3,1,1)
    plot(tmax, fmax);
    ylabel('Fundamental freq (Hz)');
    xlabel('Time (s)');

    subplot(3,1,2)
    plot(tmax, amax);
    ylabel('Fundamental amplitude');
    xlabel('Time (s)');
    
    subplot(3,1,3)
    plot(tmax, pmax);
    ylabel('Signal power');
    xlabel('Time (s)');

end

