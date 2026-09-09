clc;
clear;
close all;

[sig, Fe] = audioread('hel_fr1.wav');
sig = sig';

N          = length(sig);
L          = round(50e-3 * Fe);   % Longueur de trame (50 ms)
LW         = 2 * L;               % Longueur de la fenêtre (chevauchement 50 %)
N_trames   = floor(N / L) - 1;   
w          = sqrt(hanning(LW))';  % Fenêtre COLA

% --- 2. PARAMÈTRES DE DÉ-HÉLIUMISATION ---
alpha      = 2.5;                 % Facteur de compression (ex: 1.8 à 2.5 pour baisser les formants)
kc         = 18;                  % Coupure cepstrale pour isoler l'enveloppe sans le pitch

% Masque cepstral passe-bas (Basses quéfrences)
mask_cepstre = zeros(1, LW);
mask_cepstre(1 : kc) = 1;
mask_cepstre(end - kc + 2 : end) = 1;

% Préparation des indices pour la compression du demi-spectre
N_half   = LW / 2 + 1;            % Nombre de points de 0 à Fe/2
k_orig   = 1 : N_half;
k_source = k_orig * alpha;        % Dilatation des indices sources
k_source(k_source > N_half) = N_half; % Plafonnement à Nyquist

% --- 3. TAMPONS OLA ---
ptr = 1;
mem_synthese = zeros(1, L);
bloc_avant_fft = zeros(1, LW);
signal_filtre = zeros(1, N);

% --- 4. BOUCLE PRINCIPALE ---
for trame = 1 : N_trames
    new_frame = sig(ptr : ptr + L - 1);
    bloc_avant_fft(end/2+1 : end) = new_frame;
    
    % Fenêtrage + FFT
    xw = bloc_avant_fft .* w;
    Xf = fft(xw);
    
    % Extraction du Cepstre (Log-Magnitude)
    mag_X = abs(Xf) + 1e-6;
    log_abs_X = log(mag_X); 
    Y = fft(log_abs_X);
    
    % Isolation de l'enveloppe spectrale originale E
    Y_pb = Y .* mask_cepstre;
    E = exp(real(ifft(Y_pb))); 
    
    % Isolation de l'excitation spectrale (X / E)
    Excitation = Xf ./ E;
    
    % COMPRESSION DE L'ENVELOPPE SUR LE DEMI-SPECTRE
    E_half = E(1 : N_half);
    E_comp_half = interp1(k_orig, E_half, k_source, 'linear', 'extrap');
    
    % Reconstitution de l'enveloppe complète avec symétrie hermitienne (miroir)
    E_comp = zeros(1, LW);
    E_comp(1 : N_half) = E_comp_half;
    E_comp(N_half + 1 : end) = E_comp_half(LW/2 :-1 : 2);
    
    % Recomposition du spectre (Excitation * Enveloppe modifiée)
    Xf_mod = Excitation .* E_comp;
    
    % IFFT + Fenêtrage de synthèse + Overlap-Add
    y = real(ifft(Xf_mod));
    yw = y .* w;
    
    trame_OLA = yw(1 : end/2) + mem_synthese;
    signal_filtre(ptr : ptr + L - 1) = trame_OLA;
    
    % Mises à jour
    mem_synthese = yw(end/2+1 : end);
    bloc_avant_fft(1:end/2) = bloc_avant_fft(end/2+1:end);
    ptr = ptr + L;
end

% --- 5. AJUSTEMENT DE GAIN ET ÉCOUTE ---
signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));
sound(signal_filtre, Fe);