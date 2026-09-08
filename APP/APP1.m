% %% APP 1 Félix Boivin BOIF1302
% %% Clear and Load
% 
% clear
% [sig, Fe] = audioread('hel_fr1.wav');
% sig = sig';
% 
% %% Approche LPC
% 
% ploting = false;
% 
% clc
% N          = length(sig);
% L          = 50e-3*Fe;           % Longueur d'une trame (50 ms à Fe = 44.1 kHz)
% LW         = 2 * L;              % Longueur de la fenêtre (chevauchement de 50 %)
% N_trames   = floor(N / L) - 1;   % Nombre de trames
% w          = sqrt(hanning(LW))'; % Fenêtre d'analyse/synthèse COLA
% m          = 20;                 % Ordre LPC (m=10 à 20 recommandé)
% k          = 2.15;
% 
% ptr            = 1;
% mem_synthese   = zeros(1, L);
% mem_lpc_synth  = zeros(1, m);    % Mémoire des états pour le filtre IIR
% trame_analyse  = zeros(1, LW);
% signal_filtre  = zeros(1, N);
% 
% % boucle principale
% 
% for trame = 1 : N_trames
% 
%     % Trames
%     new_frame = sig(ptr : ptr + L - 1);
%     trame_analyse(1 : end/2) = trame_analyse(end/2+1 : end);
%     trame_analyse(end/2+1 : end) = new_frame;
% 
%     % Fenêtrage Hanning
%     xw = trame_analyse .* w;
% 
%     % Coefficients LPC
%     A_LPC = lpc(xw, m);
%     if ploting
%         %Pour grahique initiale
%         [H,w_axis] = freqz(1,A_LPC,LW/2);
%         H          = H';
%         w_axis     = w_axis';
%         Ha         = 20*log10(abs(H));
%         Xa         = 20*log10(abs(fft(xw)));
%         Xa         = Xa(1:LW/2);
%         diff_max   = 0;
% 
%             figure(1)
%             subplot(4,1,1)
%             hold off
%             plot(w_axis, Xa);
%             hold on
%             plot(w_axis,Ha - diff_max,'r')
%             title("Avant")
%             subplot(4,1,2)
%             plot(new_frame)
%             title("Signal Avant")
%     end
% 
%     % Erreur
%     err = filter(A_LPC, 1, xw);
% 
% 
%     % Manipulation des pôles
%     p = roots(A_LPC);              % Obtenir les pôles du filtre d'analyse
%     r = abs(p);                    % Rayon
%     theta = angle(p);              % Angle
% 
%     theta_mod = theta/k;
% 
%     p_mod = r .* exp(1i * theta_mod);
% 
%     A_LPC_mod= real(poly(p_mod));
% 
%     % Reconstitution signal y[n] 
%     [y, mem_lpc_synth] = filter(1, A_LPC_mod, err, mem_lpc_synth);
% 
%     % Fenêtrage Hanning 2
%     yw = y .* w;
% 
%     % OLA
%     trame_OLA = yw(1 : end/2) + mem_synthese;
%     signal_filtre(ptr : ptr + L - 1) = trame_OLA;   
%     mem_synthese = yw(end/2+1 : end);
% 
%     %Pour grahique initiale
%     if ploting
%         [H2,w_axis2] = freqz(1,A_LPC_mod,LW/2);
%         H2          = H2';
%         w_axis2     = w_axis2';
%         Ha2         = 20*log10(abs(H2));
%         Xa2         = 20*log10(abs(fft(yw)));
%         Xa2         = Xa2(1:LW/2);
%         diff_max   = 0;
% 
%             subplot(4,1,3)
%             hold off
%             plot(w_axis2, Xa2);
%             hold on
%             plot(w_axis2,Ha2 - diff_max,'r')
%             title("LPC Apres")
%             subplot(4,1,4)
%             plot(mem_synthese)
%             title("Signal apres LPC")
%     end
% 
%     ptr = ptr + L;
% 
%     if ploting
%         pause(1);
%     end
% end
% 
% signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));
% % max(abs(signal_filtre)) * max(abs(sig));
% % signal_filtre = signal_filtre / 40000;
% 
% % sound(signal_filtre, Fe);

%% APP 1 
% Félix Boivin BOIF1302
% Mathieu Desautels DESM1210

%% Clear and Load
clear
clc
[sig, Fe] = audioread('hel_fr1.wav');
sig = sig';

%% ==================== Paramètres de contrôle ====================
ploting   = false;   % Affichage des graphiques avant/après (FFT + signal)
use_LPC   = true;    % Approche 1 : LPC (filtre adaptatif à prédiction linéaire)
use_FFT   = false;   % Approche 2 : décomposition fréquentielle (FFT), sur enveloppe cepstrale
                      % Note : pour comparer les 2 approches (comme demandé dans le
                      % guide), n'active qu'un seul bloc à la fois.

N          = length(sig);
L          = 2*round(50e-3*Fe/2); % Longueur d'une trame (50 ms), forcée entière et paire
LW         = 2 * L;               % Longueur de la fenêtre (chevauchement 50 %)
N_trames   = floor(N / L) - 1;    % Nombre de trames
w          = sqrt(hanning(LW))';  % Fenêtre d'analyse/synthèse COLA

% Facteur de compression d'enveloppe spectrale (celui causé par l'hélium),
% typiquement entre 2 et 3 selon la profondeur (voir énoncé, section 4).
% On veut RETROUVER une enveloppe comprimée de ce facteur, SANS déplacer
% la position des harmoniques (donc sans changer le pitch).
comp_env = 2.5;

%% ------------- Paramètres propres au bloc LPC -------------
m = 20;   % Ordre LPC (m = 10 à 20 recommandé)

%% ------------- Paramètres propres au bloc FFT -------------
lifter_order = 30;  % Ordre du "lifter" cepstral séparant l'enveloppe (basse
                     % "quéfrence") de la structure fine/harmoniques (haute
                     % "quéfrence"). Doit rester plus petit que la période de
                     % pitch en échantillons (sinon les harmoniques se
                     % retrouvent mélangées dans l'enveloppe).

%% ==================== Initialisation ====================
ptr            = 1;
mem_synthese   = zeros(1, L);
mem_lpc_synth  = zeros(1, m);     % Mémoire des états pour le filtre IIR
trame_analyse  = zeros(1, LW);
signal_filtre  = zeros(1, N);

%% ==================== Boucle principale ====================
for trame = 1 : N_trames

    % Trames
    new_frame = sig(ptr : ptr + L - 1);
    trame_analyse(1 : end/2)     = trame_analyse(end/2+1 : end);
    trame_analyse(end/2+1 : end) = new_frame;

    % Fenêtrage Hanning (une seule fois, avant tout traitement)
    xw = trame_analyse .* w;

    %% ---------------- APPROCHE 1 : LPC ----------------
    if use_LPC

        A_LPC = lpc(xw, m);

        if ploting
            [H, w_axis] = freqz(1, A_LPC, LW/2);
            H = H'; w_axis = w_axis';
            Ha = 20*log10(abs(H));
            Xa = 20*log10(abs(fft(xw)));
            Xa = Xa(1:LW/2);

            figure(1)
            subplot(4,1,1); hold off; plot(w_axis, Xa); hold on; plot(w_axis, Ha, 'r');
            title("FFT Avant")
            subplot(4,1,2); plot(new_frame); title("Signal Avant")
        end

        % Erreur (résidu d'excitation) = signal blanchi de son enveloppe
        err = filter(A_LPC, 1, xw);

        % Manipulation des pôles pour comprimer l'enveloppe spectrale
        p     = roots(A_LPC);   % Pôles du filtre d'analyse
        r     = abs(p);         % Rayon (garde la même "netteté" de résonance)
        theta = angle(p);       % Angle (position en fréquence du formant)

        % Compression de l'enveloppe : diviser l'angle des pôles par comp_env
        % ramène chaque formant à f_orig/comp_env, ce qui comprime l'enveloppe
        % en fréquence SANS changer la structure fine (le résidu "err" -
        % où se trouvent les harmoniques - n'est pas touché).
        theta_mod = theta / comp_env;
        p_mod     = r .* exp(1i * theta_mod);
        A_LPC_mod = real(poly(p_mod));

        % Reconstitution du signal y[n] avec l'enveloppe comprimée
        [y, mem_lpc_synth] = filter(1, A_LPC_mod, err, mem_lpc_synth);

        if ploting
            [H2, w_axis2] = freqz(1, A_LPC_mod, LW/2);
            H2 = H2'; w_axis2 = w_axis2';
            Ha2 = 20*log10(abs(H2));
            Xa2 = 20*log10(abs(fft(y)));
            Xa2 = Xa2(1:LW/2);

            subplot(4,1,3); hold off; plot(w_axis2, Xa2); hold on; plot(w_axis2, Ha2, 'r');
            title("FFT Après LPC")
        end

    else
        y = xw;   % Bloc LPC désactivé : signal inchangé
    end

    %% ---------------- APPROCHE 2 : FFT (enveloppe cepstrale) ----------------
    if use_FFT

        Xf    = fft(y);
        logX  = log(abs(Xf) + eps);     % Log-magnitude du spectre

        % Séparation enveloppe / structure fine par liftrage cepstral
        cep         = real(ifft(logX));
        cep_env     = cep;
        cep_env(lifter_order+1 : end-lifter_order+1) = 0;  % ne garde que la basse quéfrence
        env_log     = real(fft(cep_env));                   % enveloppe (log-magnitude lissée)
        fine_log    = logX - env_log;                        % structure fine (contient les harmoniques)

        % Axe de fréquence équivalent à celui de freqz (0 à pi), pour être
        % directement comparable au graphique du bloc LPC
        w_axis_fft = (0:LW/2-1) * (2*pi/LW);

        if ploting
            Xa = 20*log10(abs(Xf));
            Xa = Xa(1:LW/2);
            Ha = (20/log(10)) * env_log(1:LW/2);

            figure(1)
            subplot(4,1,1);
            hold off; 
            plot(w_axis_fft, Xa); 
            hold on; 
            plot(w_axis_fft, Ha, 'r');
            title("FFT Avant")
            subplot(4,1,2); 
            plot(new_frame); 
            title("Signal Avant")
        end

        % Compression de l'enveloppe en fréquence par comp_env, en préservant
        % la symétrie conjuguée (signal réel) : on ne travaille que sur la
        % demi-bande [0, Nyquist], puis on reconstruit l'autre moitié en miroir.
        half         = LW/2;
        freq_idx     = 0:half;
        env_log_half = env_log(1:half+1);

        % env_compresse(f) = env_recue(f * comp_env)  -> ramène les formants
        % étirés par l'hélium vers leur position d'origine
        env_log_half_c = interp1(freq_idx, env_log_half, freq_idx*comp_env, 'linear', 'extrap');

        env_log_c = [env_log_half_c, fliplr(env_log_half_c(2:end-1))]; % symétrie conjuguée

        % Recombinaison : nouvelle enveloppe + structure fine d'origine (harmoniques intactes)
        logX_new = env_log_c + fine_log;
        Xf_new   = exp(logX_new) .* exp(1i * angle(Xf));

        y = real(ifft(Xf_new));

        if ploting
            Xa2 = 20*log10(abs(fft(y)));
            Xa2 = Xa2(1:LW/2);
            Ha2 = (20/log(10)) * env_log_c(1:LW/2);

            subplot(4,1,3); 
            hold off; 
            plot(w_axis_fft, Xa2); 
            hold on; 
            plot(w_axis_fft, Ha2, 'r');
            title("FFT Après")
        end
    end

    % Fenêtrage Hanning final (une seule fois, après tout traitement)
    yw = y .* w;

    % OLA
    trame_OLA = yw(1 : end/2) + mem_synthese;
    signal_filtre(ptr : ptr + L - 1) = trame_OLA;
    mem_synthese = yw(end/2+1 : end);

    if ploting
        figure(1)
        subplot(4,1,4); plot(mem_synthese); title("Signal après")
    end

    ptr = ptr + L;

    if ploting
        pause(1);
    end
end

signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));

sound(signal_filtre, Fe);