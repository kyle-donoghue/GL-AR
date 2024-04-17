clc;clear;close all hidden;
%% Importing Paths
addpath(genpath("../graph_creation/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../progressBar/"));
addpath(genpath("../comparison_tests/"));

%% Global Parameters
global_params.N = 20; % node count
global_params.l = 1e3; % signal length
global_params.noise = 1e-6; %noise power
global_params.trials = 20;
global_params.Fs = 100;
global_params.delay = 20;
global_params.minSep = 100;
global_params.maxSep = 500;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 10; % AR approx order
donoghue_params.gamma = 1.5e-5;%.04;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.interval_length = 100;
donoghue_params.threshold_avg = 0.03;
donoghue_params.threshold = .015;

%% Statistical Storage
    
%% Multi-Trial
gammas = logspace(-4,4,40);
p_orders = 5:5:45;
thresholds = [ logspace(-4,1.6,50)];
thresholds_avg = 0:.005:.5;
j_max = length(thresholds);
l_max = length(thresholds_avg);
[Fmeasure, precision, recall, NMI, num_of_edges] = deal(zeros(length(thresholds),length(thresholds_avg),length(gammas)));
count = 1;
for k = 1:length(gammas)
    for i = 1:length(p_orders)
        toc
        tic
        count/(length(gammas)*length(p_orders))*100
        count = count+1;
        donoghue_params.gamma = gammas(k);
        [f, p, r, nm, ne] = deal(zeros(length(thresholds),length(thresholds_avg),global_params.trials));
        
        parfor t = 1:global_params.trials
        
            %% Creating The Graph
            % [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
            [A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
            % [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);
            
            %% Donoghue Normalization
            A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
            L_0 = adjacency_to_laplacian(A_0);
            A_1 = to_directed(A_0);

            %% Simulating Signals
            % freqs = generate_frequencies2(global_params.N,2,25);
            freqs = 10.2:.2:14;
            X_noisy = varm_pulse_signal(A_1,global_params,freqs);
            

            %% Donoghue Graph Learning with Intervals
            L_donoghue_occ = graph_learning_AR_occ(X_noisy,donoghue_params);
            for j = 1:j_max
                for l = 1:l_max
                    A_donoghue_occ_t = laplacian_to_adjacency(L_donoghue_occ)>thresholds(j); % binarize each interval's A
                    A_donoghue_occ_avg = mean(A_donoghue_occ_t,3); % average binarized A's
                    A_donoghue_occ_avg(A_donoghue_occ_avg<thresholds_avg(l)) = 0; % threshold the averaged result
                    L_donoghue_occ_avg = adjacency_to_laplacian(A_donoghue_occ_avg); % final laplacian

                    %% Graph Accuracy
                    [p(j,l,t),r(j,l,t),f(j,l,t),nm(j,l,t),ne(j,l,t)] = graph_learning_perf_eval(L_0,L_donoghue_occ_avg);
                end
            end

       end
        
        precision(:,:,k,i) = mean(p,3);
        recall(:,:,k,i) = mean(r,3);
        Fmeasure(:,:,k,i) = mean(f,3);
        NMI(:,:,k,i) = mean(nm,3);
        num_of_edges(:,:,k,i) = mean(ne,3);
    end
end

%% Outputs
[X_axis,Y_axis] = meshgrid(thresholds,thresholds_avg);
% save("thresholds_energy_coeff_proport_occ_n20.mat");