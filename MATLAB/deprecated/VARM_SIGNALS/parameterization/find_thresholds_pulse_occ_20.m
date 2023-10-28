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
global_params.noise = .005; %noise power
global_params.trials = 20;
global_params.Fs = 100;
global_params.delay = 10;
global_params.minSep = 250;
global_params.maxSep = 500;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 50; % AR approx order
donoghue_params.gamma = .05;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.intervals = 20;
donoghue_params.interval_length = floor(global_params.l/donoghue_params.intervals);
donoghue_params.threshold_avg = 0;

%% Statistical Storage

%% Multi-Trial
thresholds_norm = .01:.02:.15;
p_orders = 5:5:20;
[Fmeasure, precision, recall, NMI, num_of_edges] = deal(zeros(length(thresholds_norm),length(p_orders),2));
for k = 1:length(thresholds_norm)
    for j = 1:length(p_orders)
        ((k-1)*length(p_orders)+j)/length(thresholds_norm)/length(p_orders)*100
        donoghue_params.threshold = thresholds_norm(k);
        donoghue_params.P = p_orders(j);
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
            freqs = generate_frequencies2(global_params.N,2,25);
            X_noisy = varm_pulse_signal(A_0,global_params,freqs);
            
            %% Donoghue Graph Learning with Intervals
            L_donoghue_occ = graph_learning_AR_occ(X_noisy,donoghue_params);
            A_donoghue_occ_t = laplacian_to_adjacency(L_donoghue_occ)>donoghue_params.threshold; % binarize each interval's A
            A_donoghue_occ_avg = mean(A_donoghue_occ_t,3); % average binarized A's
            A_donoghue_occ_avg(A_donoghue_occ_avg<donoghue_params.threshold_avg) = 0; % threshold the averaged result
            L_donoghue_occ_avg = adjacency_to_laplacian(A_donoghue_occ_avg); % final laplacian
            
            %% Graph Accuracy
            [p(t),r(t),f(t),nm(t),ne(t)] = graph_learning_perf_eval(L_0,L_donoghue_occ_avg);
            
        end
        
        precision(k,j,:) = [mean(p) std(p)];
        recall(k,j,:) = [mean(r) std(r)];
        Fmeasure(k,j,:) = [mean(f) std(f)];
        NMI(k,j,:) = [mean(nm) std(nm)];
        num_of_edges(k,j,:) = [mean(ne) std(ne)];
    end
end

%% Outputs
[X,Y] = meshgrid(thresholds_norm,p_orders);
surf(X,Y,Fmeasure(:,:,1)')
% save("thresholds_occ_n20.mat");