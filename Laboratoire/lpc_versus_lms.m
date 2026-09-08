%  LPC_versus_LMS.m
%
%  Petit script qui compare les coeffcients de prédiction obtenus
%  par l'approche LPC et par l'approche LMS (gradient stochastique)
%
%  Auteur : Roch Lefebvre
%  Date   :  24 aout 2026
%
%-------------------------------------------------------------------

[sig, Fe] = audioread('C:\Roch\Code_Octave\S7CI\APP1\Audio\voix_homme_32_kHz.wav');
sig = sig';

debut = 26500;
L = 2048;
x = sig(debut : debut + L - 1);
w = hamming(L)';

m = 2;   % ordre des filtres prédicteurs

%  Calcul des coefficients avec l'approche LPC

xw = x .* w;
A_LPC = lpc(xw, m);

%   Calcul des coefficiehts de prédiction par l'approche LMS (gradient)

alpha = 0.5;
A_LMS =   randn(1,m);

vect_erreur_LMS = zeros(1,L);
vect_erreur_LMS(1:m) = x(1:m);

memoire_A_LMS = zeros(L,2);

for ech = m+1 : length(x)
     pred = A_LMS * x(ech-1 : -1 : ech-m)';
     e    = x(ech) - pred;
     vect_erreur_LMS(ech) = e;
     A_LMS    =  A_LMS + 2 * alpha * e * x(ech-1 : -1 : ech-m);
end

err_LPC = filter(A_LPC,1,x);
%  err_LMS = filter([1 -A_LMS],1,x);   % juste si on veut filter un fois avec le foltre LMS "final"

figure(1)
hold off
plot(x)
hold on
plot(err_LPC,'r')
plot(vect_erreur_LMS,'g')
