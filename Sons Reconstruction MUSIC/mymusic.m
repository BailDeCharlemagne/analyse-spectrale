%%% Algorithme de Music pour la d\'etermination des param\`etres Music
%%%
%%% entr\'ees :
%%%   - xx : signal
%%%   - pp : ordre du mod\`ele (choisi de mani\`ere ind\'ependante)
%%%   - MM : nombre de coefficients de corr\'elation pris en compte (MM>=pp)
%%%   - fe : fr\'equence d'\'echantillonnage
%%%
%%% sorties :
%%%   - ff     : fr\'equences auxquelles la dsp a \'et\'e calcul\'ee
%%%   - mydsp  : la dsp elle-m\^eme
%%%
%%% exemples : clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1));mymusic(xx,2,10,fe);
%%%            clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1))+cos(2*pi*(f0+26)/fe*[1:tsig]+2*pi*rand(1,1));mymusic(xx,4,300,fe);
%%%
%%% S. Rossignol -- 2012

%%% Modified by Y. GREMBER -- 2016


%%% utilisation en script :
%clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1))+cos(2*pi*(f0+26)/fe*[1:tsig]+2*pi*rand(1,1));pp=4;MM=400;
%clear;rand('seed',100*sum(clock));fe=32000;f0=440;tsig=1280;xx=cos(2*pi*f0/fe*[1:tsig]+2*pi*rand(1,1))+1e-2*randn(1,tsig);pp=2;MM=10;


function [ff, mydsp] = mymusic(xx, pp, MM, fe, df, fmin, fmax)


if (MM<=pp)
  fprintf(1, 'Il faut absolument MM>pp !!!\n');
  return;
end;

MM1=MM-1;
res=xx;
xx = xx-mean(xx);


%%% corr\'elations
acf = xcorr(xx, MM1, 'biased');
lMM=length(acf); % = 2MM1 + 1
rrr1 = acf(MM1+1:lMM); % longueur MM1 + 1 %%% Penser à détransposer
for ii=1:MM1
  rrr1 = [rrr1 acf(MM1+1-ii:lMM-ii)]; %%% Penser à détransposer
end;
% rrr1 = rrr1'; % rrr1 matrice carrée de de dimension MM1 + 1
% rrr1 contient sur la ligne ii les coeffs Rxx(0-ii)..Rxx(MM1-ii)


%%% m\'ethode directe pour trouver toutes les valeurs propres
[v, lambda] = eig(rrr1);
lamb = diag(lambda);
[vl,pl] = sort(abs(lamb),'descend');


%%% densit\'e spectrale de puissance
df=0.1; %%% la dsp est calcul\'ee tous les df Hz
ff=fmin:df:fmax;

mydsp=zeros(length(ff),1);
deni=zeros(MM,MM);
for ii=pp+1:MM
  deni = deni + v(:,pl(ii))*conj(v(:,pl(ii)))';
end;
for ii=1:length(ff)
  ee = cos(2*pi*ff(ii)*[0:MM1]/fe);
  den = conj(ee)*deni*ee';
  mydsp(ii) = abs(1/den);
end;


%%% on enl\`eve \'eventuellement une composante non nulle
mydsp = mydsp-min(mydsp);
%mydsp = mydsp/max(mydsp); %%% si on fait \c{c}a, c'est norme 1


%%% figures
% figure(1);
% clf;
% grid on;
% hold on;
% plot(ff,mydsp,'linewidth',2);
% xlabel('frequency (in Hz)','fontsize',20);
% ylabel('magnitude','fontsize',20);
% xlim([300 400]);
% hold off;
% % 
% % % % % % % % if(plot != 0)
%  figure(2);
%  clf;
%  grid on;
%  hold on;
%  plot(ff,mydsp,'linewidth',2);
%  xlabel('frequency (in Hz)','fontsize',20);
%  ylabel('magnitude','fontsize',20);
%  hold off;
%  drawnow;
% % % % % % % % end if