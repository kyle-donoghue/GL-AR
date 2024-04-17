% This file serves as a general guide on how to:
% 1) Import established signal set and groundtruth graph
% 2) Set up GL-AR
% 3) Apply GL-AR to the signals and recover the groundtruth graph structure
clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
addpath(genpath('../Formal_Generation/new_temperature'));
load('temperature_dataset.mat')

%% Signal/Graph Import
[L_dist,A_dist] = create_temperature_groundtruth(DIST,100,15); % using DISTance matrix, calculate adjacency with sigma, alpha parameters
density = .01; % assume edge density to be 1%
[y_noisy,~,sigma] = signals.z_score(detrended(:,1:2:end)); % z-score detrended temperature data

%% GL-AR Parameter Creation
signal_params = signals.create_empty(size(y_noisy)); % create a blank signal_params struct based off size of data
AR_params = GL.create_default_params(signal_params); % default params for GL-AR with a given signal size (P=40 by default)
AR_params.P = 64; % change to P=64 (found manually)
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
[p,r,f,nmi,n_e] = graphs.performance(L_dist,L_t); % calculate precision, recall, f-score, nmi, number of edges
fprintf("F-score: %.2f%%\n", f*100)
