clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
signal_params.N = 8;
signal_params.raw = 1;
%%
L_0 = graphs.create(signal_params,'er',0.2);
A_0 = graphs.to_adjacency(L_0);
%%
N = signal_params.N;
M = 100;

w = randn(N,M);

order = 1;
c = [0 .25];

x = w;

for k = 1+order:M
    for i = 1:order
        total_A = zeros(N);
        for j = 0:i
            total_A = total_A + c(j+1)*A_0^j;
        end
        x(:,k) = x(:,k) + total_A*x(:,k-i);
    end
end

plot(x')