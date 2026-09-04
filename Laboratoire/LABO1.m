
% Import des sons
fprintf('Eagles\n')
[eagles_48k, Fs1] = audioread('eagles_48k.wav');
fprintf('yellow_48k\n')
[yellow_48k, Fs7] = audioread('yellow_48k.wav');

%Figures de départ
figure()
subplot(2,1,1)
plot(eagles_48k)
title('Eagles 48k')
subplot(2,1,2)
plot(yellow_48k)
title('Yellow 48k')

help hanning()

