clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_approximation/"));
%% Parameters
N = 4;
M=100;
num_active = 5;
P = 26;
gamma = 200;
%% Creating The Graph
A = zeros(N);
A(1,3) = 1;
A(2,4) = .3;
A = 1/2*(A+A');

%% Trace Normalization
L_0 = adjacency_to_laplacian(A);
L_0 = L_0*N/trace(L_0);

%% Signal Creation
load("../signal_creation/testmat.mat")
figure
plot(out')
V = out(:,1:1:end);
V = V(:,1:2000)';
% V = out(:,1:20:end);
% V = V(:,1:100);
%% Filter
% original 250Hz cutoff with 2000Hz Fs
% pulse takes ~10-20ms
% out(:,1:2000) ~ 50ms
% 2000*(1/.05) is new Fs = 40e3
% create filter with Fs 40e3 and 250Hz cutoff:
% 250Hz -> 2*pi*250/40e3 = .0125*pi
figure
subplot(1,2,1)
signalPlot(V)
filt = firpm(200,[0 .0125 .025 1],[1 1 0 0]);
V = filter(filt,1,V);
subplot(1,2,2)
signalPlot(V)
%% Z score
V = (V-mean(V))./std(V);
%% Approx
B = approxAR(V,P);
% B = B-mean(B);
%% Graph Learning
Q = create_Q_matrix(N,1);
c = create_c_vec(N,B')*gamma;
A = create_constraint_matrix(N);
b = create_constraint_vec(N,N);
options = optimset('Display', 'off');
[phi,f] = quadprog(Q,c,A(N+2:end,:),b(N+2:end),A(1:N+1,:),b(1:N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
l = create_dup_matrix(N)*phi;
L = convert_to_matrix(l);
% L(abs(L)<.05) = 0;

%% P Order
aics = [];
for P = 2:2:50
    aics(end+1) = ar(V(:,1)',P,'ls').Report.Fit.AIC;
end
figure
plot(2:2:50,aics)
%% Metrics
figure
subplot(1,2,1)
signalPlot(V)
subplot(1,2,2)
signalPlot(B)
figure
vcompare(L_0,L)
% [precision, recall, Fmeasure, NMI, num_of_edges] = graph_learning_perf_eval(L_0,L)