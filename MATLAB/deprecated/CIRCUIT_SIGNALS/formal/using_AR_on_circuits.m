clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
%%
N = 20;
M = 1e3;
P = 20;
num_active = 20;
noise = 1.8e-4;
shunt = 8;
divide_resist = 60;
Fs = 100;
T = 1/Fs:1/Fs:M/Fs;

gamma = 2.6;
threshold = .085;

%%
[A_0,XCoords, YCoords] = construct_graph(N,'er',0.2);
% [A,XCoords, YCoords] = construct_graph(N,'pa',1);
L_0 = adjacency_to_laplacian(A_0);
L_0 = L_0*N/divide_resist/trace(L_0);
L_0 = L_0+shunt*eye(N);
I = zeros(N,M);
active = randperm(N,num_active);
freqs = generate_frequencies2(num_active,2,45);
I(active,:) = sin(2*pi*freqs*T);
V = pinv(L_0)*I+noise*randn(N,M);
V = V';
V = (V-mean(V))./std(V);
V = V';
B = approxAR(V',P);
Q = create_Q_matrix(N,1);
c = create_c_vec(N,B')*gamma;
A = create_constraint_matrix(N);
b = create_constraint_vec(N,N);
options = optimset('Display', 'off');
[phi,f] = quadprog(Q,c,A(N+2:end,:),b(N+2:end),A(1:N+1,:),b(1:N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
l = create_dup_matrix(N)*phi;
L = convert_to_matrix(l);
L(abs(L)<threshold) = 0;
A_new = laplacian_to_adjacency(L);
%% Metrics
[p,r,f,nmi,num_edges] = graph_learning_perf_eval(L_0,L)
vcompare(A_0,A_new)