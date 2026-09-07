%   visualisaton_enveloppe_LPC.m
%
%   Script Matlab qui affiche l'enveloppe d'un filtre LPC
%   (enveloppe du filtre de synthèse 1/A(z)) et la superpose
%   au spectre d'amplitude de la trame du signal analysé.
%
%   On boucle sur plusieurs trames consécutives du signal
%   pour montrer que l'enveloppe du filtre de synthèse
%   s'adapte (i.e. modélise) l'enveloppe du signal
%   (à un gain près)
%
%---------------------------------------

[sig, Fe] = audioread('C:\Roch\Code_Octave\S7CI\APP1\Audio\voix_homme_32_kHz.wav');
sig = sig';

N = length(sig);     % Nbre d'échantillons du sihnal complet
L = 2*320;            % Longueur de la fenêtre d'analyse
w = hamming(L)';     % Fenêtre d'analyse
N_trames = floor(N/L)-1;    % nombre de trames de longueur L pleines

m = 10;              % ordre du filtre LPC

premiere_trame = 50;
derniere_trame = 80;

ptr = (premiere_trame-1) * L + 1;

for trame = premiere_trame : derniere_trame
    x  = sig(ptr : ptr + L - 1);
    xw = x .* w;
    a  = lpc(xw,m);
    [H,w_axis] = freqz(1,a,L/2);
    H          = H';
    w_axis     = w_axis';
    Ha         = 20*log10(abs(H));
    Xa         = 20*log10(abs(fft(xw)));
    Xa         = Xa(1:L/2);

    diff_max   = 0 % max(Ha) - max(Xa);   %  mettre à 0 pour voir la diff entre les spectres

        figure(1)
        hold off
        plot(w_axis,Xa);
        hold on
        plot(w_axis,Ha - diff_max,'r')
        figure(2)
        plot(x)

    ptr = ptr + L;
    pause(1)
end

