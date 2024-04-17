clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
rng("default")
%% Parameters
N = 8;
M=1e4;
num_active = 1;
Fs = 100;
T = 1/Fs:1/Fs:M/Fs;

%% Creating The Graph
% [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
[A,XCoords, YCoords] = construct_graph(N,'er',0.3);
% [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);
figure
plot(graph(A))
%% Trace Normalization
L_0 = adjacency_to_laplacian(A);
L_0 = L_0*N/50/trace(L_0);

%% Signal Creation
I = zeros(N,M);
active = randperm(N,num_active);
% for i = 1:num_active
%     L_0(active(i),active(i)) = L_0(active(i),active(i))+.1;
% end
% L_0(2,2) = L_0(2,2)+.1;
L_0 = L_0+.1*eye(N);
freqs = generate_frequencies2(num_active,2,10);
I(active,:) = cos(2*pi*freqs*T);
% freqs = generate_frequencies2(N,2,10);
% I = sin(2*pi*freqs*T);

figure
subplot(1,3,1)
signalPlot(I')
subplot(1,3,2)
V = pinv(L_0)*I+1e-2*randn(N,M);
signalPlot(V')
V = V';
V = (V-mean(V))./std(V);
V = V';
subplot(1,3,3)
signalPlot(V')
%% Graph Learning
gamma = 50;
B = approxAR(V',20);
Q = create_Q_matrix(N,1);
c = create_c_vec(N,B')*gamma;
A = create_constraint_matrix(N);
b = create_constraint_vec(N,N);
options = optimset('Display', 'off');
[phi,f] = quadprog(Q,c,A(N+2:end,:),b(N+2:end),A(1:N+1,:),b(1:N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
l = create_dup_matrix(N)*phi;
L = convert_to_matrix(l);
L(abs(L)<.075) = 0;

%% Metrics

figure
signalPlot(B)
vcompare(L_0,L)
[precision, recall, Fmeasure, NMI, num_of_edges] = graph_learning_perf_eval(L_0,L)

figure
aics = [];
for p=2:100
    aics(end+1) = ar(V(3,:),p,'ls').Report.Fit.AIC;
end
plot(2:100,aics)