%% APP 1 Félix Boivin BOIF1302
%% Clear and Load

clear
[sig, Fe] = audioread('hel_fr1.wav');
sig = sig';
%% Approche LPC

clc
N          = length(sig);
L          = 50e-3*Fe;                % Longueur d'une trame (50 ms à Fe = 44.1 kHz)
LW         = 2 * L;              % Longueur de la fenêtre (chevauchement de 50 %)
N_trames   = floor(N / L) - 1;    % Nombre de trames
w          = sqrt(hanning(LW))'; % Fenêtre d'analyse/synthèse COLA
m          = 20;                 % Ordre LPC (m=10 à 20 recommandé)
k          = 2.15;

ptr            = 1;
mem_synthese   = zeros(1, L);
mem_lpc_synth  = zeros(1, m);    % Mémoire des états pour le filtre IIR
trame_analyse  = zeros(1, LW);
signal_filtre  = zeros(1, N);

% boucle principale

for trame = 1 : N_trames

    % Trames
    new_frame = sig(ptr : ptr + L - 1);
    trame_analyse(1 : end/2) = trame_analyse(end/2+1 : end);
    trame_analyse(end/2+1 : end) = new_frame;
    
    % Fenêtrage Hanning
    xw = trame_analyse .* w;
    
    % Coefficients LPC
    A_LPC = lpc(xw, m);
    
    % Erreur
    err = filter(A_LPC, 1, xw);

    % Manipulation des pôles
    p = roots(A_LPC);              % Obtenir les pôles du filtre d'analyse
    r = abs(p);                    % Rayon
    theta = angle(p);              % Angle

    theta_mod = theta/k;

    p_mod = r .* exp(1i * theta_mod);

    A_LPC_mod= real(poly(p_mod));
    
    % Reconstitution signal y[n] 
    [y, mem_lpc_synth] = filter(1, A_LPC_mod, err, mem_lpc_synth);
    
    % Fenêtrage Hanning 2
    yw = y .* w;
    
    % OLA
    trame_OLA = yw(1 : end/2) + mem_synthese;
    signal_filtre(ptr : ptr + L - 1) = trame_OLA;   
    mem_synthese = yw(end/2+1 : end);
    
    ptr = ptr + L;
end

% signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));
% max(abs(signal_filtre)) * max(abs(sig));
signal_filtre = signal_filtre / 40000;

sound(signal_filtre, Fe);