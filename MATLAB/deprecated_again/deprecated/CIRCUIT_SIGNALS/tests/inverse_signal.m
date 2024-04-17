clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
%%
N = 8;
M = 100;
%% Creating The Graph
% [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
[A,XCoords, YCoords] = construct_graph(N,'er',0.2);
% [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);

%% Trace Normalization
L_0 = adjacency_to_laplacian(A);
L_0 = L_0*N/trace(L_0);

%%
[V,D] = eig(L_0);
B = V*mvnrnd(zeros(N,1),pinv(D),M)';

%%
b = B(1,:);
h = freqz(1,b,4096,'whole');
x = ifft(h);
figure
plot(abs(h))
figure
plot(x)