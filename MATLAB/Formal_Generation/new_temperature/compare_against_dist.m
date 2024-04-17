clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('temperature_dataset.mat')

%% Find distance adjacency

%%
% good_indices = [125 24 25 126 23 124 115 70 130 7 87 49 119 138 21 97 20 93 57 145 63 94 14 16  88 103 73 82 39 72 141 96 70 17 149 117 59 145 52 130 86 135 78 64 34 35 104 67 102 46];
% good_indices = sort(good_indices,'ascend');
% DIST = DIST(good_indices,good_indices);
%%
sigma = 100;
alpha = 15;
A_dist = zeros(length(DIST));
for m = 1:length(DIST)
    for n = 1:length(DIST)
        if m==n
            continue;
        end
        A_dist(m,n) = exp(-(DIST(m,n)/sigma)^2);
        denom1 = 0;
        denom2 = 0;
        [sort_row,ind1] = sort(DIST(m,:),'ascend');
        [sort_col,ind2] = sort(DIST(:,n),'ascend');
        for i = 1:alpha
            denom1 = denom1+exp(-(DIST(m,ind1(i))/sigma)^2);
            denom2 = denom2+exp(-(DIST(ind2(i),n)/sigma)^2);
        end
        A_dist(m,n) = A_dist(m,n)/sqrt((denom1*denom2));
    end
end
L_dist = graphs.to_laplacian(A_dist);
L_dist = GL.threshold(L_dist,1e-2);
A_dist = graphs.to_adjacency(L_dist);
graphs.plot(A_dist,'Adist');
dist_density = graphs.density(L_dist)
%% Problem Parameters
% density = .3;
density = dist_density;
[y_noisy,~,sigma] = signals.z_score(detrended(:,1:2:end));
signal_params = signals.create_empty(size(y_noisy));
N = size(y_noisy,1);
%% Solve GL-AR
AR_params = GL.create_default_params(signal_params);
AR_params.gamma = GL.get_gamma(density);
max_edges = graphs.max_edges(signal_params.N);
num_edges = ceil(density*max_edges);
pred_edges = max(1,num_edges-1);
AR_params.P = 64;
L = GL.AR_mean(y_noisy,AR_params);
weights = graphs.get_weights(L);
t = GL.get_threshold(weights,pred_edges,density);
L_t = GL.threshold(L,t);
[p,r,f,n,~] = graphs.performance(L_dist,L_t)

graphs.vcompare(A_dist,graphs.to_adjacency(L_t),'Adist','A');

%% Solve GL-SigRep
dong_params = signal_params;
dong_params.l = dong_params.interval_length;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.lambda = 10.^1;
betas = logspace(-4,4,24);
thresholds = .001:.001:1;
% thresholds = .0001;
t_max = length(thresholds);
Fd = zeros(length(betas),length(thresholds));
pd = zeros(length(betas),length(thresholds));
rd = zeros(length(betas),length(thresholds));
nd = zeros(length(betas),length(thresholds));
parfor i = 1:length(betas)
    i
    
    % dong_params.beta = betas(i);
    [Ld,~,~] = GL.dong(y_noisy(:,1:signal_params.interval_length),dong_params,betas(i));
    for t = 1:t_max
        Ld = GL.threshold(Ld,thresholds(t));
        [pd(i,t),rd(i,t),Fd(i,t),nd(i,t),~] = graphs.performance(L_dist,Ld);
    end
end
[maxi,ind] = max(Fd,[],'all');
[row,col] = ind2sub(size(Fd),ind);
pd(row,col)
rd(row,col)
Fd(row,col)
nd(row,col)