clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_approximation/"));
addpath(genpath("../signal_creation/"));
%% Parameters
N = 20;
M=100;
num_active = 5;
trials = 20;
P = 50;
Fs = 100;

gammas = logspace(-8,8,40);
thresholds = 0:1e-4:.25;
t_max = length(thresholds);
F = zeros(length(gammas),length(thresholds));
T = 1/Fs:1/Fs:M/Fs;

for k = 1:length(gammas)
    fmeasures = zeros(trials,length(thresholds));
    toc
    tic
    k/length(gammas)*100
    gamma = gammas(k);
    for i = 1:trials
        [A,XCoords, YCoords] = construct_graph(N,'er',0.2);
        L_0 = adjacency_to_laplacian(A);
        L_0 = L_0*N/trace(L_0);
        
        I = zeros(N,M);
        active = randperm(N,num_active)
        % I(active,:) = randn(num_active,M);
        freqs = generate_frequencies2(num_active,2,10);
        I(active,:) = sin(2*pi*freqs*T);
        V = pinv(L_0)*I+1e-3*randn(N,M);
        
        Q = create_Q_matrix(N,1);
        B = approxAR(V',P);
        c = create_c_vec(N,B')*gamma;
        A = create_constraint_matrix(N);
        b = create_constraint_vec(N,N);
        options = optimset('Display', 'off');
        [phi,f] = quadprog(Q,c,A(N+2:end,:),b(N+2:end),A(1:N+1,:),b(1:N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
        l = create_dup_matrix(N)*phi;
        L = convert_to_matrix(l);
        
        for t = 1:t_max
            L_tmp = L;
            L_tmp(abs(L_tmp)<thresholds(t)) = 0;
            [~,~,fmeasures(i,t),~,~] = graph_learning_perf_eval(L_0,L_tmp);
        end
    end
    F(k,:) = mean(fmeasures,1);
end
%% Metrics
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
s =  surf(X_axis,Y_axis,F');
s.EdgeColor = 'none';
