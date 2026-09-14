%% Question 1
clc 
clear all
close all

N = 2000;
m = 4;

p1 = 0.95*exp(0.1*pi*j);
p2 = conj(p1);

p3 = 0.98*exp(0.45*pi*j);
p4 = conj(p3);

a1 = conv([1 -p1], [1 -p2]);
a2 = conv([1 -p3], [1 -p4]);

a = conv(a1, a2);

e = randn(1,N);

x = filter(1, a, e);

w = hamming(N)';

xw = x.*w;

A = lpc(xw,m);

R0 = sum(xw.*xw);
R1 = sum(xw(1:end-1).*xw(2:end));
R2 = sum(xw(1:end-2).*xw(3:end));
R3 = sum(xw(1:end-3).*xw(4:end));
R4 = sum(xw(1:end-4).*xw(5:end));

R = [R0 R1 R2 R3; R1 R0 R1 R2; R2 R1 R0 R1; R3 R2 R1 R0];
r = [R1; R2; R3; R4];

A_sol = inv(R)*r;

A_sol
A

%% Question 2
clc
clear all
close all
N = 10000;
e = randn(1, N);
h = [zeros(1,20) ones(1,5) zeros(1,20) -ones(1,5)];

y = filter(h, 1, e);

%% Question 3
clc
clear all
close all

Fe = 16000;
N = floor(1*Fe);
n = 0: N-1

x = sin(2*pi*60*n/Fe)+sin(2*pi*1000*n/Fe);

L = 256;
LW = L*2;
N_trames = floor(Fe/N)-1;
w = sqrt(hanning(LW))';

ptr = 1;
mem_synthese = zeros(1, L);
bloc_avant_fft = zeros(1, LW);
bloc_apres_ifft = zeros(1,LW);
signal_filtre = zeros(1, N);

for trames = 1 : N_trames
   new_frame = x(ptr:ptr+L-1); 

   bloc_avant_fft(end/2+1:end) = new_frame;
   xw = bloc_avant_fft .* w;
end

%% Question 1
clc
clear all
close all

N = 2000;
m = 4;

p1 = 0.95*exp(0.1*pi*j);
p2 = conj(p1);

p3 = 0.98*exp(0.45*pi*j);
p4 = conj(p3);

a1 = real(conv([1 -p1], [1 -p2]));
a2 = real(conv([1 -p3], [1 -p4]));
a = conv(a1, a2);

e = randn(1, N);

x = filter(1, a, e);

w = hamming(N)';

xw = x .* w;

A = lpc(xw, m);

R0 = sum(xw.*xw);
R1 = sum(xw(1:end-1).*xw(2:end));
R2 = sum(xw(1:end-2).*xw(3:end));
R3 = sum(xw(1:end-3).*xw(4:end));
R4 = sum(xw(1:end-4).*xw(5:end));

R = [R0 R1 R2 R3; R1 R0 R1 R2; R2 R1 R0 R1; R3 R2 R1 R0];
r = [R1;R2;R3;R4];

A_sol = inv(R)*r;

%% Question 2
clc
clear all
close all

N = 10000;
h = [zeros(1,20) ones(1,5) zeros(1,20) -ones(1,5)];

x = randn(1, N);
m = 60;
A_LMS = randn(1,m);
Alpha = 0.001;

y = filter(h,1,x);

for ech = m+1 : N
    pred = A_LMS * x(ech-1:-1:ech-m)';
    e = y(ech) - pred;
    A_LMS = A_LMS + 2*Alpha*e*x(ech-1:-1:ech-m);
end

figure(1)
plot(A_LMS)

%% Question 2
clc
clear all
close all

N = 10000;
x = randn(1, N);

h = [zeros(1,20) ones(1,5) zeros(1,20) -ones(1,5)];

Alpha = 10/N;
m = 60;

y = filter(h, 1, x);
A_LMS = randn(1,m);

for ech = m+1:N
    pred = A_LMS * x(ech-1:-1:ech-m)';
    e = y(ech)-pred;
    A_LMS = A_LMS + 2*Alpha*x(ech-1:-1:ech-m)*e;
end

figure(1)
plot(A_LMS)


%% Question 1
clc
clear all
close all

N = 2000;
m = 4;

p1 = 0.95*exp(0.1*j*pi);
p2 = conj(p1);

p3 = 0.98*exp(0.45*pi*j);
p4 = conj(p3);

a1 = real(conv([1 -p1], [1 -p2]));
a2 = real(conv([1 -p3], [1 -p4]));
a = conv(a1,a2);

e = randn(1,N);
x = filter(1,a,e);

w = hamming(N)';

xw = x .* w;

A = lpc(xw, m);

