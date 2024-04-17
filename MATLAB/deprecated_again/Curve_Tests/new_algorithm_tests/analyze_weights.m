clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('weights_n20.mat');
load('weights_n10.mat');
%%
figure;

W = W_ER_20;
D = D_ER_20;
T = T_ER_20;
graph_param = er_param;
graph_type = 'er';

% W = W_ER_10;
% D = D_ER_10;
% T = T_ER_10;
% graph_param = er_param;
% graph_type = 'er';

% W = W_RND_10;
% D = D_RND_10;
% T = T_RND_10;
% graph_param = rnd_param;
% graph_type = 'gaussian';
%%
for i = 1:length(graph_param)
    density = graphs.get_sparsity(graph_type,graph_param(i));
    num_edges = ceil(density*signal_params.N*(signal_params.N-1)/2);
    for j = 1:trials
        subplot(1,trials,j)
        histogram(-W(:,j,i));
        xline(T(j,i))
        sorted_weights = sort(-W(:,j,i),'descend');
        cutoff_weight = sorted_weights(num_edges);
        xline(cutoff_weight,'Color','magenta')
        title(density)
    end
    sorted_weights = sort(-W(:,3,i),'descend');
    cutoff_weight = sorted_weights(num_edges)
    T(3,i)
end