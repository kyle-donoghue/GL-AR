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

sparsities = S_ER_20;
sparse_F = F_ER_20;
graph_type = 'er';
graph_param = er_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
cliff_delay = zeros(length(sparsities),1);
edge_spread = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    

    [a,b] = GL.get_cliff_curve(sparsities(i));
    ridge = .25;
  
    edge_spread(i) = graphs.get_edge_spread(graph_type,graph_param(i),N);

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    
    cliff_delay(i) = (log(a)+b*z_g-log(z_t-ridge))/b;
end
plot(edge_spread,real(cliff_delay));hold on;



sparsities = S_RND_20;
sparse_F = F_RND_20;
graph_type = 'gaussian';
graph_param = rnd_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
cliff_delay = zeros(length(sparsities),1);
edge_spread = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    

    [a,b] = GL.get_cliff_curve(sparsities(i));
    ridge = .25;
  
    edge_spread(i) = graphs.get_edge_spread(graph_type,graph_param(i),N);

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    
    cliff_delay(i) = (log(a)+b*z_g-log(z_t-ridge))/b;
end
plot(edge_spread,real(cliff_delay));



sparsities = S_BA_20;
sparse_F = F_BA_20;
graph_type = 'pa';
graph_param = ba_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
cliff_delay = zeros(length(sparsities),1);
edge_spread = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    

    [a,b] = GL.get_cliff_curve(sparsities(i));
    ridge = .25;
  
    edge_spread(i) = graphs.get_edge_spread(graph_type,graph_param(i),N);

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    
    cliff_delay(i) = (log(a)+b*z_g-log(z_t-ridge))/b;
end
plot(edge_spread,real(cliff_delay));

figure;
N = 10;

sparsities = S_ER_10;
sparse_F = F_ER_10;
graph_type = 'er';
graph_param = er_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
cliff_delay = zeros(length(sparsities),1);
edge_spread = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    

    [a,b] = GL.get_cliff_curve(sparsities(i));
    ridge = .25;
  
    edge_spread(i) = graphs.get_edge_spread(graph_type,graph_param(i),N);

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    
    cliff_delay(i) = (log(a)+b*z_g-log(z_t-ridge))/b;
end
plot(edge_spread,real(cliff_delay));hold on;



sparsities = S_RND_10;
sparse_F = F_RND_10;
graph_type = 'gaussian';
graph_param = rnd_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
cliff_delay = zeros(length(sparsities),1);
edge_spread = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    

    [a,b] = GL.get_cliff_curve(sparsities(i));
    ridge = .25;
  
    edge_spread(i) = graphs.get_edge_spread(graph_type,graph_param(i),N);

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    
    cliff_delay(i) = (log(a)+b*z_g-log(z_t-ridge))/b;
end
plot(edge_spread,real(cliff_delay));



sparsities = S_BA_10;
sparse_F = F_BA_10;
graph_type = 'pa';
graph_param = ba_param;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
cliff_delay = zeros(length(sparsities),1);
edge_spread = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    

    [a,b] = GL.get_cliff_curve(sparsities(i));
    ridge = .25;
  
    edge_spread(i) = graphs.get_edge_spread(graph_type,graph_param(i),N);

    [z,ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    z_g = log10(gammas(row));
    z_t = thresholds(col)*N;
    
    cliff_delay(i) = (log(a)+b*z_g-log(z_t-ridge))/b;
end
plot(edge_spread,real(cliff_delay));