clc
clear all
close all

% Import des sons
[Sig1, Fe] = audioread('eagles_48k.wav');
[Sig2, Fe] = audioread('yellow_48k.wav');

%Figures de départ
figure()
subplot(2,1,1)
plot(Sig1)
title('Eagles 48k')
subplot(2,1,2)
plot(Sig2)
title('Yellow 48k')
close

%Information
N_Sig1 = length(Sig1);
N_Sig2 = length(Sig2);
L = 1024; %Dit dans la question
LW = 2 * L; %Longueur de la fenetre
N_trames_Sig1 = floor(N_Sig1/L) - 1;
N_trames_Sig2 = floor(N_Sig2/L) - 1;
w = sqrt(hanning(LW));

mask_LP = [1 ones(1,LW/4-1) 0 zeros(1, LW/2) ones(1, LW/4-1)];
mask_HP = 1 - mask_LP;
    ind_min = floor((300 / Fe) * LW);
    ind_max = floor((3400 / Fe) * LW);
    n_ones = ind_max - ind_min + 1;
mask_BP = zeros(1, LW);
mask_BP(ind_min: ind_max) = ones(1, n_ones);
mask_BP(LW - ind_max: LW - ind_min) = ones(1, n_ones);

ptr = 1;
mem_synthese_Sig1 = zeros(1, L);
mem_synthese_Sig2 = zeros(1, L);
bloc_avant_fft_Sig1 = zeros(1, LW);
bloc_avant_fft_Sig2 = zeros(1, LW);
bloc_apres_ifft_Sig1 = zeros(1, LW);
bloc_apres_ifft_Sig2 = zeros(1, LW);
signal_filtre_Sig1 = zeros(1, N_Sig1);
signal_filtre_Sig2 = zeros(1, N_Sig2);

%Pour Sig2
for trame = 1 : N_trames_Sig2
        new_frame = Sig2(ptr : ptr + L-1);
        bloc_avant_fft_Sig2(end/2+1 : end) = new_frame;
        xw = bloc_avant_fft_Sig2 .* w;
        Xf = fft(xw);

        Xf_mod = Xf .* mask_BP;

        y = real(ifft(Xf_mod));
        yw = y .* w;
        trame_OLA = yw(1 : end/2) + mem_synthese_Sig2;
        signal_filtre(ptr : ptr + L -1) = trame_OLA;

        mem_synthese_Sig2 = yw(end/2+1 : end);
        bloc_avant_fft_Sig2(1:end/2) = bloc_avant_fft_Sig2(end/2+1:end);
        ptr = ptr + L;
end

%Pour Sig1
for trame = 1 : N_trames_Sig1
        new_frame = Sig1(ptr : ptr + L-1);
        bloc_avant_fft_Sig1(end/2+1 : end) = new_frame;
        xw = bloc_avant_fft_Sig1 .* w;
        Xf = fft(xw);

        Xf_mod = Xf .* mask_BP;

        y = real(ifft(Xf_mod));
        yw = y .* w;
        trame_OLA = yw(1 : end/2) + mem_synthese_Sig1;
        signal_filtre(ptr : ptr + L -1) = trame_OLA;

        mem_synthese_Sig1 = yw(end/2+1 : end);
        bloc_avant_fft_Sig1(1:end/2) = bloc_avant_fft_Sig1(end/2+1:end);
        ptr = ptr + L;
end
