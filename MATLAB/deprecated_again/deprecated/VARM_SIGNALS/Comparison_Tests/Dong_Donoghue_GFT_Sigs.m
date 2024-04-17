clc;clear;close all hidden;
%% Importing Paths
addpath(genpath("../graph_creation/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
addpath(genpath("../GL_algorithms/"));
%% Global Parameters
global_params.N = 20; % node count
global_params.l = 100; % signal length
global_params.noise = .5; %noise power
global_params.trials = 1;

%% Dong Parameters
dong_params = global_params;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.beta = 10.^(-0.2);
dong_params.lambda = 10.^1;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 30; % AR approx order
donoghue_params.gamma = .15;
donoghue_params.threshold = .065;

%% Multi-Trial
for t = 1:global_params.trials
    t
    %% Creating The Graph
    % [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
    [A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
    % [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);
    
    %% Dong Normalization
    L = adjacency_to_laplacian(A);
    L_0 = L/trace(L)*global_params.N; % make trace of L_0=N

    %% Simulating Signals
    X_noisy = dong_gft_signal(L_0,global_params);
    
    %% Signal Approximation
    B = approxAR(X_noisy',donoghue_params.P);
    
    %% Donoghue Graph Learning
    L_donoghue = graph_learning_AR(B',donoghue_params);
    L_donoghue(abs(L_donoghue)<donoghue_params.threshold) = 0;
    %% Dong Graph Learning
    [L_dong,Y,~] = graph_learning_gaussian(X_noisy,dong_params);
    L_dong(abs(L_dong)<10^(-4))=0; % prune insignificant weights
    
    %% Graph Accuracy
    [precision(t,1),recall(t,1),Fmeasure(t,1),NMI(t,1),num_of_edges(t,1)] = graph_learning_perf_eval(L_0,L_dong);
    [precision(t,2),recall(t,2),Fmeasure(t,2),NMI(t,2),num_of_edges(t,2)] = graph_learning_perf_eval(L_0,L_donoghue);
end
%% Outputs
mean(Fmeasure,1)