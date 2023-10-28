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
global_params.noise = .0; %noise power
global_params.trials = 20;
global_params.Fs = 100;
global_params.freqs = generate_frequencies2(global_params.N,2,25);
global_params.delay = 25;
global_params.minSep = 50;
global_params.maxSep = 100;

%% Dong Parameters
dong_params = global_params;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.beta = 10.^(-0.2);
dong_params.lambda = 10.^1;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 30; % AR approx order
donoghue_params.gamma = .05;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.threshold = .08;
donoghue_params.intervals = 5;
donoghue_params.interval_length = floor(global_params.l/donoghue_params.intervals);

%% Statistical Storage
[Fmeasure, precision, recall, NMI, num_of_edges] = deal(zeros(global_params.trials,3));

%% Multi-Trial
for 
    parfor t = 1:global_params.trials
    
        %% Creating The Graph
        % [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
        [A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
        % [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);
        
        %% Donoghue Normalization
        A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
        L_0 = adjacency_to_laplacian(A_0);
    
        %% Simulating Signals
        freqs = generate_frequencies2(global_params.N,2,25);
        % freqs = zeros(20,1);
        X_noisy = varm_pulse_signal(A_0,global_params,freqs);
        
        %% Signal Approximation
        B = approxAR(X_noisy',donoghue_params.P);
    
        %% Donoghue Graph Learning
        L_donoghue_0 = graph_learning_AR(B',donoghue_params);
        L_donoghue = L_donoghue_0;
        L_donoghue(abs(L_donoghue)<donoghue_params.threshold) = 0;
    
        %% Donoghue Graph Learning with Intervals
        L_donoghue_occ = graph_learning_AR_occ(X_noisy,donoghue_params);
        A_donoghue_occ_t = laplacian_to_adjacency(L_donoghue_occ)>donoghue_params.threshold;
        L_donoghue_occ_avg = adjacency_to_laplacian(mean(A_donoghue_occ_t,3));
        
        %% Dong Graph Learning
        [L_dong,Y,~] = graph_learning_gaussian(X_noisy,dong_params);
        L_dong(abs(L_dong)<10^(-4))=0; % prune insignificant weights
        
        %% Graph Accuracy
        [p,r,f,nm,ne] = deal(zeros(1,3));
        [p(1),r(1),f(1),nm(1),ne(1)] = graph_learning_perf_eval(L_0,L_dong);
        [p(2),r(2),f(2),nm(2),ne(2)] = graph_learning_perf_eval(L_0,L_donoghue);
        [p(3),r(3),f(3),nm(3),ne(3)] = graph_learning_perf_eval(L_0,L_donoghue_occ_avg);
        precision(t,:) = p;
        recall(t,:) = r;
        Fmeasure(t,:) = f;
        NMI(t,:) = nm;
        num_of_edges(t,:) = ne;
        
    end
end
%% Outputs
figure(1)
histogram(Fmeasure(:,1))
figure(2)
histogram(Fmeasure(:,2))
figure(3)
histogram(Fmeasure(:,3))
mean(Fmeasure,1)