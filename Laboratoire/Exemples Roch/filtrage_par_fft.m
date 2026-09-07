%  filtrage_par_FFT.m
%
%  Démonstration de filtrage par transformée (ici, la FFT).
%  Approche par trames, avec fenêtrage complémentaire à 50 %.
%
%  Attention : Il faut fenêtrer AVANT la FFT et APRÈS la FFT inverse.
%              Il faut donc que le produit des fenêtres (produit point à point)
%              donne une fenêtre complémentaire pour recouvrement à 50 %
%              dans l'overlap-and-add.
%
%              Ici, on va prendre la fenêtre de Hanning, et mettre chaque élément
%              à la puissance 1/2 (racine carrée) pour formée la fenêtre qui
%              sera utilisée avant la FFT et après la FFT inverse.
%
% Auteur : Roch Lefebvre


[sig, Fe] = audioread('C:\Roch\Code_Octave\S7CI\APP1\Audio\eagles_48k.wav');
sig = sig';

N          = length(sig);
L          = 1024;      %  longueur d'une trame (20 ms à Fe = 48 kHz)
LW         = 2 * L;        % longueur de la fenêtre (2 trames, avec overla de 50 %)
N_trames   = floor(N / L) - 1;  % nombre de trames à parcourir
w          = sqrt(hanning(LW))';  % fenêtre d'analyse et de synthèse

%  Les masques pour le filtrage fréquentiel

  mask_LP = [1 ones(1,LW/4-1) 0 zeros(1,LW/2) ones(1,LW/4-1)];
  mask_HP = 1 - mask_LP;
     ind_min = floor((300 / Fe) * LW);
     ind_max = floor((3400 / Fe) * LW);
     n_ones = ind_max - ind_min + 1;
  mask_BP = zeros(1,LW);
  mask_BP(ind_min:ind_max) = ones(1,n_ones);
  mask_BP(LW-ind_max:LW-ind_min) = ones(1,n_ones);


ptr = 1;
mem_synthese = zeros(1,L);
bloc_avant_fft = zeros(1,LW);
bloc_apres_ifft = zeros(1,LW);
signal_filtre = zeros(1,N);

% boucle principale

for trame = 1 : N_trames
  new_frame = sig(ptr : ptr + L - 1);
  bloc_avant_fft(end/2+1 : end) = new_frame;
  xw = bloc_avant_fft .* w;
  Xf = fft(xw);

  Xf_mod = Xf .* mask_HP;     % on modifie le spectre ici

  y = real(ifft(Xf_mod));
  yw = y .* w;
  trame_OLA = yw(1 : end/2) + mem_synthese;
  signal_filtre(ptr : ptr + L - 1) = trame_OLA;   % on va écrire cette nouvelle trame dans le tableau de sortie global

  mem_synthese = yw(end/2+1 : end);
  bloc_avant_fft(1:end/2) = bloc_avant_fft(end/2+1:end);
  ptr = ptr + L;
end
