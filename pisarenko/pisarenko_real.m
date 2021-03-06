## Cette fonction renvoie uniquement la partie réelle et 
## un vecteur de pics de fréquence à partir de la méthode Pisarenko 
#       xx = signal vecteur 
#       order = ordre de la méthode pisarenko 
#       freq_sampling = freq échantillonnage du signal 
#       option =  0 si on veut juste les vpropres et la dsp 
#                 1 si on veut afficher le spectre de raies 
#                 2 si on veut l'enveloppe interpolée 
#       plot_or_not = 'plot' si on veut afficher la courbe

function [ff_real, mydsp_real] = pisarenko_real(xx, order, freq_sampling, option, plot_or_not)
  
  [ff_pisa, mydsp_pisa] = mypisarenko(xx,order,freq_sampling, option); % donne densité de puissance spectrale
  ff_real = ff_pisa(floor(length(ff_pisa)/2)+1:end); 
  mydsp_real = mydsp_pisa(floor(length(mydsp_pisa)/2)+1:end);
 
  if (logical(plot_or_not) == 1)
    figure; 
    plot(ff_real, mydsp_real);
    str=sprintf('Plot with parameter %d', option);
    title(str);
  end;
end; 
  