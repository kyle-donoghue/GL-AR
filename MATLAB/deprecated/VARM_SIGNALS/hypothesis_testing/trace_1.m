clc;clear;close all hidden;
%% Importing Paths
addpath(genpath("../graph_creation/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../progressBar/"));
%% Global Parameters
global_params.N = 20; % node count
global_params.l = 1e3; % signal length
global_params.noise = .05; %noise power
global_params.trials = 20;
global_params.Fs = 250;
global_params.freqs = generate_frequencies2(global_params.N,5,100);
global_params.delay = 50;

[A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
L = adjacency_to_laplacian(A);
L = L/trace(L)*10; % set trace to 1;

A_0 = laplacian_to_adjacency(L);

X_noisy = varm_signal(A_0,global_params,global_params.freqs);
A_0(:,1)
plot(abs(fft(X_noisy(1,:))))
max(eig(A_0))
