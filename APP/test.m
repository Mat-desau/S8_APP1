%% test comp
clear
clc

c = mod(10,3);

n = 100;
k = 2.7;

vect = rand(1,n);
new_vect = [];

pad = []

for i = 1 : n

    j = mod(i,k);
    disp(j);

    if j == 0
        new_vect =[new_vect, vect(i)];
    end
    
end

pad = zeros(1,n-length(new_vect));
new_vect = [new_vect, pad]