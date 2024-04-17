clc;clear;close all hidden;
load('large_maximum_data_results_n20_n10.mat');
load('large_maximum_data_results_n40.mat');
addpath(genpath('../../GL_classes/'));
%%
% sparsities = S_RND_20;
% sparse_F = F_RND_20(:,2:end,:);
% N = 20;
% graph_type = 'gaussian';
% graph_param = rnd_param;
% thresholds = 0:.0025:2;
% gammas = logspace(-2,4,20);
%%
sparsities = S_ER_20;
sparse_F = F_ER_20(:,2:end,:);
N = 20;
graph_type = 'er';
graph_param = rnd_param;
thresholds = 0:.0025:2;
gammas = logspace(-2,4,20);
%%
calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
f_t = zeros(length(sparsities),1);
f_g = zeros(length(sparsities),1);
f_z = zeros(length(sparsities),1);
for i = 1:(length(sparsities))
    [f_z(i), ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    f_t(i) = thresholds(col);
    f_g(i) = gammas(row);

end
f_t = f_t*N;
f_g = log10(f_g);
%%
% f_g = f_g(1:7);
% f_g(7) = 0;
% f_t = f_t(1:7);
%%
figure;hold on;
plot(f_t,f_g);
%%
f_ti = f_t(end):.1:f_t(1);
% f_gi = interp1(f_t,f_g,f_ti);
% figure;hold on;
% plot(f_ti,f_gi)
%%
a = 1.4756;
b = -3.5826;
c = 0.8592;
d = 0.0340;
y_hat = a+b*exp(-c*(f_ti-d)).*cos(-c*(f_ti-d));
plot(f_ti,y_hat)