% This file serves as a general guide on how to:
% 1) Establish synthetic signal/graph parameters
% 2) Create a synthetic signal set using the Graph Tensor Method (GTM)
% 3) Apply GL-AR to the signals and recover the groundtruth graph structure
clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));

%% Signal/Graph Parameter Creation
N = 20; % 20 vertices
graph_type = 'er'; % Erdos-Reyni graph
graph_param = .2; % ER graph with 20% edge density
signal_params = signals.create_default(N,graph_type); % default parameters defined in signals class

%% Graph Creation
[L_0,~,A_d] = graphs.create(signal_params,graph_param); % create ER graph and return adjacency/laplacian
G = graphs.createGraphTensor(signal_params,A_d); % create graph tensor from groundtruth adjacency

%% Signal Creation
y_noisy = signals.generateFilteredRectPulse(signal_params,G); % create noisy signal set with GTM
y_noisy = signals.z_score(y_noisy); % z-score signal set

%% GL-AR Parameter Creation
AR_params = GL.create_default_params(signal_params); % default params for GL-AR with a given signal size (P=40 by default)
density = graphs.get_density(signal_params.graph_type,graph_param,signal_params.N); % assume the edge density of a given random graph type
AR_params.gamma = GL.get_gamma(density); % use automated parameter selection algorithm to choose gamma based on edge density

%% GL-AR
L = GL.AR_mean(y_noisy,AR_params); % estimate autoregressive matrices for each signal interval, average the autoregression matrices, and solve the GL-AR cost function

max_edges = graphs.max_edges(signal_params.N); % maximum amount of edges for a given graph size
num_edges = ceil(density*max_edges); % expected edges based on assumed edge density
pred_edges = max(1,num_edges-1); % take away one from num_edges

weights = graphs.get_weights(L); % form edge weight vector
t = GL.get_threshold(weights,pred_edges,density); % calculate corresponding threshold using automated parameter selection algorithm

L_t = GL.threshold(L,t); % threshold raw recovered graph laplacian

%% Evaluation
[p,r,f,nmi,n_e] = graphs.performance(L_0,L_t); % calculate precision, recall, f-score, nmi, number of edges
fprintf("F-score: %.2f%%\n", f*100)
