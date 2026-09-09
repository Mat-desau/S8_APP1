%% APP 1 Félix Boivin BOIF1302
%% Clear and Load

clear
[sig, Fe] = audioread('hel_fr1.wav');
sig = sig';

%% Approche LPC

clc
N          = length(sig);
L          = 50e-3*Fe;           % Longueur d'une trame (50 ms à Fe = 44.1 kHz)
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
        title("FFT Avant")
        subplot(4,1,2)
        plot(new_frame)
        title("Signal Avant")
    
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
        title("FFT Apres")
        subplot(4,1,4)
        plot(mem_synthese)
        title("Signal apres")
    
    ptr = ptr + L;
    pause(1)
end

signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));
% max(abs(signal_filtre)) * max(abs(sig));
% signal_filtre = signal_filtre / 40000;

sound(signal_filtre, Fe);

%% Approche FFT

clc
close all

N          = length(sig);
L          = 50e-3*Fe;           % Longueur d'une trame (50 ms à Fe = 44.1 kHz)
LW         = 2 * L;              % Longueur de la fenêtre (chevauchement de 50 %)
N_trames   = floor(N / L) - 1;   % Nombre de trames
w          = sqrt(hanning(LW))'; % Fenêtre d'analyse/synthèse COLA
m          = 20;                 % Ordre LPC (m=10 à 20 recommandé)
k          = 3;

% Passe-bas
mask_PB = [1 ones(1,floor(LW/4)) 0 zeros(1,LW/2) ones(1,floor(LW/4)-1)];
% plot(mask_PB)

ptr = 1;
mem_synthese = zeros(1, L);
bloc_avant_fft = zeros(1, LW);
signal_filtre = zeros(1, N);

% Boucle principale
for trame = 1 : N_trames
    new_frame = sig(ptr : ptr + L - 1);
    bloc_avant_fft(end/2+1 : end) = new_frame;
    
    % fenetre hanning
    xw = bloc_avant_fft .* w;

    % premiere fft
    Xf = fft(xw);

    log_abs_X = log(abs(Xf) + 1e-6); % + 1e-6 évite log(0)
    Y = fft(log_abs_X);
    
    % On ne garde que les premiers et derniers indices (quefrences lentes)
    kc = 18; 
    mask_cepstre = zeros(1, LW);
    mask_cepstre(1 : kc) = 1;
    mask_cepstre(end - kc + 2 : end) = 1;
    
    Y_pb = Y .* mask_cepstre;
    
    % 5. Retour au domaine fréquentiel pour obtenir l'enveloppe E
    log_E = real(ifft(Y_pb));
    E = exp(log_E); % Exponentielle pour annuler le log de l'étape 3

    plot(E(1:end/2))

    E_temp = [];
    pad = [];

    for i = 1 : LW
        j = mod(i,k);
        disp(j);
        if j == 0
            E_temp = [E_temp, E(i)];
        end
    end

    E_temp_half = E_temp(1:floor(length(E_temp));


    %pad = zeros(1,LW-length(E_comp));
    %E_comp = [E_comp, pad]

    figure;
    plot(E_comp)
    
    % % 5. IFFT + Fenêtrage de synthèse
    % y = real(ifft(Xf_mod));
    % yw = y .* w;
    % 
    % % 6. Overlap-And-Add (OLA)
    % trame_OLA = yw(1 : end/2) + mem_synthese;
    % signal_filtre(ptr : ptr + L - 1) = trame_OLA;
    % 
    % % Mises à jour des tampons
    % mem_synthese = yw(end/2+1 : end);
    % bloc_avant_fft(1:end/2) = bloc_avant_fft(end/2+1:end);
    % ptr = ptr + L;
    pause(5)
end

% Normalisation d'amplitude finale pour éviter la saturation
signal_filtre = signal_filtre / max(abs(signal_filtre)) * max(abs(sig));

sound(signal_filtre, Fe);