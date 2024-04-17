clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('pollution_data.mat')
%%
sigma = 50;
alpha = 10;
A_dist = zeros(length(dist_adj));
for m = 1:length(dist_adj)
    for n = 1:length(dist_adj)
        if m==n
            continue;
        end
        A_dist(m,n) = exp(-(dist_adj(m,n)/sigma)^2);
        denom1 = 0;
        denom2 = 0;
        [sort_row,ind1] = sort(dist_adj(m,:),'ascend');
        [sort_col,ind2] = sort(dist_adj(:,n),'ascend');
        for i = 1:alpha
            denom1 = denom1+exp(-(dist_adj(m,ind1(i))/sigma)^2);
            denom2 = denom2+exp(-(dist_adj(ind2(i),n)/sigma)^2);
        end
        A_dist(m,n) = A_dist(m,n)/sqrt((denom1*denom2));
    end
end
L_dist = graphs.to_laplacian(A_dist);
L_dist = GL.threshold(L_dist,1e-2);
A_dist = graphs.to_adjacency(L_dist);
graphs.plot(A_dist,'Adist');
dist_density = graphs.density(L_dist)
%%
dataset = pollution;
dataset = fillmissing(dataset,'linear',2,'SamplePoints',1:366);
density = graphs.density(L_dist);


%%
signal_params.N = size(dataset,1);
signal_params.interval_length = size(dataset,2);
signal_params.intervals = 1;
signal_params.M = size(dataset,2);
AR_params = GL.create_default_params(signal_params);
y = signals.z_score(dataset);
%%
% [fits,P] = GL.sweep_fits(pollution);
% plot(P,fits)
%%
AR_params.gamma = GL.get_gamma(density);
AR_params.P = 65;
max_edges = graphs.max_edges(signal_params.N);
num_edges = ceil(density*max_edges);
pred_edges = max(1,num_edges-1);
L = GL.AR_mean(y,AR_params);

weights = graphs.get_weights(L);
t = GL.get_threshold(weights,pred_edges,density);

L_tmp = GL.threshold(L,t);
[f1,f2,f3,f4,f5] = graphs.performance(L_dist,L_tmp);
graphs.vcompare(L_dist,L,'L_0','Ltmp')
f3
%%
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
    
    % dong_params.beta = betas(i);
    [Ld,~,~] = GL.dong(y,dong_params,betas(i));
    for t = 1:t_max
        Ld = GL.threshold(Ld,thresholds(t));
        [~,~,Fd(i,t),~,~] = graphs.performance(L_dist,Ld);
    end
end
max(Fd,[],"all")