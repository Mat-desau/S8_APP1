
fprintf('----- UPLOAD SONS -----\n')
fprintf('hel_fr1\n')
[hel_fr1, Fs2] = audioread('hel_fr1.wav');
fprintf('hel_fr2\n')
[hel_fr2, Fs3] = audioread('hel_fr2.wav');
fprintf('hel_fr4\n')
[hel_fr4, Fs4] = audioread('hel_fr4.wav');
fprintf('parole\n')
[parole, Fs5] = audioread('parole.wav');
fprintf('parole_2\n')
[parole_2, Fs6] = audioread('parole_2.wav');


fprintf('----- GRAPHIQUES SONS -----\n')
figure()
%Hel_fr1
subplot(5,1,1)
plot(hel_fr1)
title('Hel Fr1')
%Hel_fr2
subplot(5,1,2)
plot(hel_fr2)
title('Hel Fr2')
%Hel_fr4
subplot(5,1,3)
plot(hel_fr4)
title('Hel Fr4')
%Parole
subplot(5,1,4)
plot(parole)
title('Parole')
%Parole2
subplot(5,1,5)
plot(parole_2)
title('Parole 2')








