clc;clear;close all hidden;
rng("default")
%% Importing Paths
addpath(genpath("../graph_creation/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../progressBar/"));
addpath(genpath("../phd_code/"));
addpath(genpath("../comparison_tests/"));

%% Global Parameters
global_params.N = 20; % node count
global_params.l = 1e3; % signal length
global_params.noise = 1e-6; %noise power
global_params.trials = 20;
global_params.Fs = 100;
global_params.delay = 10;
global_params.minSep = 100;
global_params.maxSep = 500;

%% Donoghue Parameters
donoghue_params = global_params;
donoghue_params.P = 10; % AR approx order
donoghue_params.gamma = 2.8e-2;%1.5e-5;%.04;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.threshold_avg = 0;%0.03;
donoghue_params.threshold = 3.7;%.015;

%changing for this test only
donoghue_params.interval_length = 100;%global_params.l;
donoghue_params.l = 100;

%% Creating The Graph
% [A,XCoords, YCoords] = construct_graph(global_params.N,'gaussian',0.75,0.5);
[A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
% [A,XCoords, YCoords] = construct_graph(global_params.N,'pa',1);

%% Donoghue Normalization
A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
L_0 = adjacency_to_laplacian(A_0);
A_1 = to_directed(A_0);

%% Simulating Signals
freqs = generate_frequencies2(global_params.N,2,25);
%freqs = 10.2:.2:14;
X_noisy = varm_pulse_signal(A_1,global_params,freqs);
X_noisy = X_noisy(:,50:149);
%% Raw Donoghue
L1 = graph_learning_AR_occ(X_noisy,donoghue_params);
A1 = laplacian_to_adjacency(L1(:,:,2));

%% Raw PhD
A2 = gl_ar(X_noisy,donoghue_params.gamma,donoghue_params.P);

%% Compare
signalPlot(X_noisy');
mse = mean((A2(:)-A1(:)).^2)
vcompare(A1,A2);
figure
imagesc(A_0)