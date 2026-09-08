clc
clear all
close all

% Import des sons
[Sig1, Fe] = audioread('eagles_48k.wav');

%Information
N = length(Sig1);
L = 1024; %Dit dans la question
N_trames = N/(L - 1);
w = hanning(L);

