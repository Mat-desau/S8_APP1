%% APP 1 Félix Boivin BOIF1302
%% Clear and Load

clear
[sig, Fe] = audioread('hel_fr1.wav');
sig = sig';

%% Approche LPC

ploting = false;

clc
N          = length(sig);
L          = 2*round(50e-3*Fe/2);           % Longueur d'une trame (50 ms à Fe = 44.1 kHz)
LW         = 2 * L;              % Longueur de la fenêtre (chevauchement de 50 %)
N_trames   = floor(N / L) - 1;   % Nombre de trames
w          = sqrt(hanning(LW))'; % Fenêtre d'analyse/synthèse COLA
m          = 20;                 % Ordre LPC (m=10 à 20 recommandé)
k          = 2.15;

ptr            = 1;
mem_synthese   = zeros(1, L);
mem_lpc_synth  = zeros(1, m);    % Mémoire des états pour le filtre IIR
trame_analyse  = zeros(1, LW);
signal_filtre  = zeros(1, N);

% =========================================================
%    LPC
% =========================================================

for trame = 1 : N_trames

    % Trames
    new_frame = sig(ptr : ptr + L - 1);
    trame_analyse(1 : end/2) = trame_analyse(end/2+1 : end);
    trame_analyse(end/2+1 : end) = new_frame;

    % Fenêtrage Hanning
    xw = trame_analyse .* w;

    % Coefficients LPC
    A_LPC = lpc(xw, m);
    if ploting
        %Pour grahique initiale
        [H,w_axis] = freqz(1,A_LPC,LW/2);
        H          = H';
        w_axis     = w_axis';
        Ha         = 20*log10(abs(H));
        Xa         = 20*log10(abs(fft(xw)));
        Xa         = Xa(1:LW/2);
        diff_max   = 0;

            figure(1)
            subplot(4,1,1)
            hold off
            plot(w_axis, Xa);
            hold on
            plot(w_axis,Ha - diff_max,'r')
            title("Avant")
            subplot(4,1,2)
            plot(new_frame)
            title("Signal Avant")
    end

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

    %Pour grahique initiale
    if ploting
        [H2,w_axis2] = freqz(1,A_LPC_mod,LW/2);
        H2          = H2';
        w_axis2     = w_axis2';
        Ha2         = 20*log10(abs(H2));
        Xa2         = 20*log10(abs(fft(yw)));
        Xa2         = Xa2(1:LW/2);
        diff_max   = 0;

            subplot(4,1,3)
            hold off
            plot(w_axis2, Xa2);
            hold on
            plot(w_axis2,Ha2 - diff_max,'r')
            title("LPC Apres")
            subplot(4,1,4)
            plot(mem_synthese)
            title("Signal apres LPC")
    end

    ptr = ptr + L;

    if ploting
        pause(1);
    end
end

signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));
% max(abs(signal_filtre)) * max(abs(sig));
% signal_filtre = signal_filtre / 40000;
fprintf("ici1")
% sound(signal_filtre, Fe);
% pause(5)

% =========================================================
%    FFT
% =========================================================
sig = sig';

ploting_cepstre = true;   % true = affiche l'enveloppe cepstrale à chaque trame
K = 10;                    % Ordre de troncature cepstrale

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
signal_filtre = zeros(1,N);

if ploting_cepstre
    figure(2)
end

% boucle principale
for trame = 1 : N_trames
    new_frame = sig(ptr : ptr + L - 1);
    bloc_avant_fft(end/2+1 : end) = new_frame;

    % Fenêtrage
    xw = bloc_avant_fft .* w;

    % FFT (une seule fois, la phase est conservée pour la resynthèse)
    Xf = fft(xw);

    % --- Visualisation optionnelle de l'enveloppe cepstrale ---
    if ploting_cepstre
        Xa = abs(Xf);
        Y2 = fft(Xa);
        Y2(K : end-K+2) = 0;   % on garde seulement les basses "quéfrences"
        E  = real(ifft(Y2));

        figure(2)
        hold on
        plot(abs(E(1 : end/2)));
        title(sprintf('Enveloppe cepstrale - Trame %d, K = %d', trame, K));
        xlabel('Échantillons'); ylabel('|E|');
        ylim([0, max(abs(E))*1.1 + eps]);
        drawnow;
        pause(0.05);
    end

    % --- Filtrage réel (masquage + resynthèse) ---
    Xf_mod = Xf .* mask_HP;     % on modifie le spectre ici

    y = real(ifft(Xf_mod));
    yw = y .* w;
    trame_OLA = yw(1 : end/2) + mem_synthese;
    signal_filtre(ptr : ptr + L - 1) = trame_OLA;
    mem_synthese = yw(end/2+1 : end);

    bloc_avant_fft(1 : end/2) = bloc_avant_fft(end/2+1 : end);
    ptr = ptr + L;
end

fprintf("ici2")
% sound(signal_filtre, Fe);

