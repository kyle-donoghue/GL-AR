clc;clear;close all hidden;
rng("default")
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
global_params.noise = 1e-4; %noise power
global_params.trials = 20;
global_params.Fs = 100;
global_params.delay = 20;
global_params.minSep = 100;
global_params.maxSep = 500;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 50; % AR approx order
donoghue_params.gamma = 500;%1.5e-5;%.04;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.interval_length = 100;% making 1d < length < 2d
donoghue_params.threshold_avg = 0;%0.03;
donoghue_params.threshold = 3.7;%.015;
    
%% Creating The Graph
% [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
[A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
% [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);

%% Donoghue Normalization
A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
L_0 = adjacency_to_laplacian(A_0);
A_1 = to_directed(A_0);
figure
matrixPlot(A_0)
figure
%% Simulating Signals
freqs = generate_frequencies2(global_params.N,2,25);
%freqs = 10.2:.2:14;
X_noisy = varm_pulse_signal(A_1,global_params,freqs);


%% Donoghue Graph Learning with Intervals
L_donoghue_occ = graph_learning_AR_occ_fft(X_noisy,donoghue_params);
A_donoghue_occ_t = laplacian_to_adjacency(L_donoghue_occ)>donoghue_params.threshold; % binarize each interval's A
A_donoghue_occ_avg = mean(A_donoghue_occ_t,3); % average binarized A's
A_donoghue_occ_avg(A_donoghue_occ_avg<donoghue_params.threshold_avg) = 0; % threshold the averaged result
L_donoghue_occ_avg = adjacency_to_laplacian(A_donoghue_occ_avg); % final laplacian

%% Graph Accuracy
[precision, recall, Fmeasure, NMI, num_of_edges] = graph_learning_perf_eval(L_0,L_donoghue_occ_avg)

vcompare(L_0,L_donoghue_occ_avg)
figure
signalPlot(X_noisy')