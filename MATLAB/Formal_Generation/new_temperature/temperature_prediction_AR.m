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
[y_noisy,~,sigma] = signals.z_score(detrended(:,1:2:end));
base = detrended(:,2:2:end);
signal_params = signals.create_empty(size(y_noisy));
N = size(y_noisy,1);
M = 2;
%% Solve GL-AR
AR_params = GL.create_default_params(signal_params);
AR_params.P = 36;

AR_params.gamma = GL.get_gamma(density);
max_edges = graphs.max_edges(signal_params.N);
num_edges = ceil(density*max_edges);
pred_edges = max(1,num_edges-1);

L = GL.AR_mean(y_noisy,AR_params);
weights = graphs.get_weights(L);
t = GL.get_threshold(weights,pred_edges,density);
L_t = GL.threshold(L,t);

A = graphs.to_adjacency(L_t);

%% Solve GL-SigRep
dong_params = signal_params;
dong_params.l = dong_params.interval_length;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.beta = 10.^(-0.2);
dong_params.lambda = 10.^1;
[Ld,~,~] = GL.dong(y_noisy(:,1:signal_params.interval_length),dong_params);
Ld = GL.threshold(Ld,10^-4);

Ad = graphs.to_adjacency(Ld);
%%
[MSE,x,c,worst] = predict_temperature(y_noisy,A,base,sigma);
[MSEd,xd,cd] = predict_temperature(y_noisy,Ad,base,sigma);
%%
MSE
MSEd
worst
subplot(1,3,1)
plot(base')
subplot(1,3,2)
plot(x')
subplot(1,3,3)
plot(xd')

%% Evaluate Graphs
graphs.vcompare(A_dist,A,'Adist','AR')
[~,~,f,~,~] = graphs.performance(graphs.to_laplacian(A_dist),graphs.to_laplacian(A))
[~,~,fd,~,~] = graphs.performance(graphs.to_laplacian(A_dist),graphs.to_laplacian(Ad))
%% functions
function [MSE,x,c,worst] = predict_temperature(X,A,Y,sigma)
    M=2;
    %% Find Matrices for c
    y_A = vec(get_x(X,M,M)-A*get_x(X,M-1,M));
    B_A = [];
    for i = 2:M
        for j = 0:i
            B_A(:,end+1) = vec(A^j*get_x(X,M-i,M));
        end
    end
    %% solve for c
    n = 3;
    gamma_2 = 1;
    cvx_begin
        variable c(n)
        minimize ( 1/2*sum(sum_square(y_A - B_A*c)) + gamma_2*norm(c,1) )
    cvx_end
    c = [0;-.5;c];
    %% Predict x_hat
    K = size(Y,2);
    N = size(Y,1);
    for k = 1:K
        x(:,k) = sigma.*randn(N,1);
        count = 1;
        for i = 1:M
            if (k-i)>0
                h = 0;
                for j=0:i
                    h = h+c(count)*A^j;
                    count = count+1;
                end
                x(:,k) = x(:,k) + h*x(:,k-i);
            end
        end
    end
    MSE = 1/(N*(K-M))*sum(norm(Y(:,M+1:K)-x(:,M+1:K))^2);
    worst = 1/(N*(K-M))*sum(norm(Y(:,M+1:K)-0*x(:,M+1:K))^2);
end
function X = get_x(x,m,M)
    X = x(:,m+1:m+size(x,2)-M);
end