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
global_params.noise = .5; %noise power
global_params.trials = 50;
global_params.Fs = 250;
global_params.freqs = generate_frequencies2(global_params.N,5,100);
global_params.delay = 50;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 30; % AR approx order
donoghue_params.gamma = .05;
donoghue_params.eig_padding = .9; % what should the max eigenvalue of A be
donoghue_params.intervals = 5;
donoghue_params.interval_length = floor(global_params.l/donoghue_params.intervals);

%% Statistical Storage

%% Multi-Trial
thresholds = .01:.01:.3;
gammas = logspace(-3,1,50);
[Fmeasure, precision, recall, NMI, num_of_edges] = deal(zeros(length(thresholds),length(gammas),2));
for k = 1:length(thresholds)
    for j = 1:length(gammas)
        ((k-1)*length(gammas)+j)/length(thresholds)/length(gammas)*100
        donoghue_params.threshold = thresholds(k);
        donoghue_params.gamma = gammas(j);
        [f, p, r, nm, ne] = deal(zeros(1,global_params.trials));
    
        parfor t = 1:global_params.trials
        
            %% Creating The Graph
            % [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
            [A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
            % [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);
            
            %% Donoghue Normalization
            A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
            L_0 = adjacency_to_laplacian(A_0);
        
            %% Simulating Signals
            X_noisy = varm_signal(A_0,global_params);
            
            %% Signal Approximation
            B = approxAR(X_noisy',donoghue_params.P);
        
            %% Donoghue Graph Learning
            L_donoghue_0 = graph_learning_AR(B',donoghue_params);
            L_donoghue = L_donoghue_0;
            L_donoghue(abs(L_donoghue)<donoghue_params.threshold) = 0;
            
            %% Graph Accuracy
            [p(t),r(t),f(t),nm(t),ne(t)] = graph_learning_perf_eval(L_0,L_donoghue);
            
        end
        
        precision(k,j,:) = [mean(p) std(p)];
        recall(k,j,:) = [mean(r) std(r)];
        Fmeasure(k,j,:) = [mean(f) std(f)];
        NMI(k,j,:) = [mean(nm) std(nm)];
        num_of_edges(k,j,:) = [mean(ne) std(ne)];
    end
end

%% Outputs
save("threshold_gamma_n20.mat");