clc;clear;close all hidden;
load('large_maximum_data_results_n15_n10_3.mat');
load('large_maximum_data_results_n20_3.mat');
addpath(genpath('../../GL_classes/'));
%%
N = 20;

% sparsities = S_RND_20;
% sparse_F = F_RND_20;
% graph_type = 'gaussian';
% graph_param = rnd_param;

% sparsities = S_ER_20;
% sparse_F = F_ER_20;
% graph_type = 'er';
% graph_param = er_param;

sparsities = S_BA_20;
sparse_F = F_BA_20;
graph_type = 'pa';
graph_param = ba_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
figure(1);hold on;
figure(2);hold on;
for i = 1:(length(sparsities))
    
    [z,~] = max(sparse_F(:,:,i),[],'all');
    
    [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds*N);
    figure(1);
    s = surf(X_axis,Y_axis,sparse_F(:,:,i)'/z);
    s.EdgeColor = "none";
    view(2)
    ylim([0,20]) 

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    scatter3(z_g,z_t,1);
    
    xlabel("Gamma")
    ylabel("Threshold")

    figure(2);
    s = surf(X_axis,Y_axis,sparse_F(:,:,i)');
    s.EdgeColor = "none";
    view(2)
    ylim([0,20])   
end