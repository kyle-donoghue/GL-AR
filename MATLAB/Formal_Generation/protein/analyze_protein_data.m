clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../../Formal_Data/protein/');
load('protein_adjacencies.mat')
%%
dataset = readmatrix('1. cd3cd28.xls')';
L_0 = graphs.to_laplacian(full_protein_adjacency);
weights = graphs.get_weights(L_0);
density = nnz(weights)/length(weights);
%%
signal_params.N = size(dataset,1);
signal_params.interval_length = floor(size(dataset,2)/4);
signal_params.intervals = 4;
signal_params.M = size(dataset,2);
AR_params = GL.create_default_params(signal_params);
y = signals.z_score(dataset);
%%
% [fits,P] = GL.sweep_fits(dataset);
% plot(P,fits)
%%
AR_params.gamma = GL.get_gamma(density);
AR_params.P = 15;
max_edges = graphs.max_edges(signal_params.N);
num_edges = ceil(density*max_edges);
pred_edges = max(1,num_edges-1);
L = GL.AR_mean(y+0*randn(size(y)),AR_params);

weights = graphs.get_weights(L);
t = GL.get_threshold(weights,pred_edges,density);

L_tmp = GL.threshold(L,t);
[f1,f2,f3,f4,f5] = graphs.performance(L_0,L_tmp);
figure;
graphs.vcompare(L_0,L,'L_0','Ltmp')
f3
%%
dong_params.N = signal_params.N;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.beta = 10.^(.7);
dong_params.threshold = 10^-1;

[L_dong,Y,~] = GL.dong(y,dong_params);
L_dong = GL.threshold(L_dong,dong_params.threshold);
[f1,f2,f3,f4,f5] = graphs.performance(L_0,L_dong);
figure;
graphs.vcompare(L_0,L_dong,'L','Ldong')
f3