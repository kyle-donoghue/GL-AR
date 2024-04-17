clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('temperature_dataset.mat')

%% Find distance adjacency
sigma = 200;
alpha = 10;
A_dist = zeros(150);
for m = 1:150
    for n = 1:150
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
dist_density = graphs.density(L_dist)
%% Problem Parameters
% density = .3;
density = dist_density;
detrended = interp1(1:365,detrended',1:.25:(365-.25),"spline")';
[y_noisy,~,sigma] = signals.z_score(detrended);
signal_params = signals.create_empty(size(y_noisy));
N = size(y_noisy,1);
%% Solve GL-AR
AR_params = GL.create_default_params(signal_params);
AR_params.P = 100;

AR_params.gamma = GL.get_gamma(density);
max_edges = graphs.max_edges(signal_params.N);
num_edges = ceil(density*max_edges);
pred_edges = max(1,num_edges-1);

L = GL.AR_mean(y_noisy,AR_params);
weights = graphs.get_weights(L);
t = GL.get_threshold(weights,pred_edges,density);
L_t = GL.threshold(L,t);
[~,~,f,~,~] = graphs.performance(L_dist,L_t)

%% Solve GL-SigRep
dong_params = signal_params;
dong_params.l = dong_params.interval_length;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.lambda = 10.^1;
betas = logspace(-4,4,40);
thresholds = .001:.001:1;
t_max = length(thresholds);
Fd = zeros(length(betas),length(thresholds));
parfor i = 1:length(betas)
    i
    
    % dong_params.beta = betas(i);
    [Ld,~,~] = GL.dong(y_noisy(:,1:signal_params.interval_length),dong_params,betas(i));
    for t = 1:t_max
        Ld = GL.threshold(Ld,thresholds(t));
        [~,~,Fd(i,t),~,~] = graphs.performance(L_dist,Ld);
    end
end
