%  function [y,ind] = quant_scal_unif(x, val_min, val_max, n_bits)
%
%  Fonction qui prend en entrée un signal numérique x
%  ainsi que les paramètres d'un quantificateur scalaire uniforme.
%
%  La fonction retourne :
%
%      - le signal quantifié y (on peut entrer tout un bloc d'échantioons dans la fonction)
%      - la série d'indices binaires ind associés à chaque niveau de reconstruction pour
%        chaque échantillon de x en entrée
%
%  Auteur : Roch Lefebvre
%

function [y, ind] = quant_scal_unif(x, val_min, val_max, n_bits)

  n_niveaux = 2 ^ n_bits;
  delta = (val_max - val_min) / n_niveaux;
  un_sur_delta = 1 / delta;

  ind = round((x - val_min) * un_sur_delta);
    ind = max(ind, 0);              % saturation "en bas"  (min indice = 0)
    ind = min(ind, n_niveaux - 1);  % saturation "en haut"  (max indice = 2^n_bits - 1)
  y   = val_min + ind * delta;

end
