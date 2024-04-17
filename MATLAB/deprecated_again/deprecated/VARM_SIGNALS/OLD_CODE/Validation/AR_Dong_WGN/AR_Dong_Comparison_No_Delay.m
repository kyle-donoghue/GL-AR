clc;clear;close all hidden;
%% AR Parameters
P = 30;
alpha = 1;
beta = .15;

%% STATISTICAL STORAGE
trials = 5e2;
f_scores = zeros(trials,1);

%% Dong Parameters
tolerance = 1e-3;

%% Global Parameters
trials = 5e2;
f_scores_AR = zeros(trials,1);
f_scores_Dong = zeros(trials,1);
precision_AR = zeros(trials,1);
precision_Dong = zeros(trials,1);
recall_AR = zeros(trials,1);
recall_Dong = zeros(trials,1);

n = 20;
noise = .25;
N = 2e3;
delay = 50;
interval = 1;
Fs = 250;
freqs = generate_frequencies2(n,5,100);
ER_prob = .2;

threshold = .065;
%% Data Set Generation
[coeffs, y, x, A] = sineER(n,P,ER_prob,delay,interval,N,Fs,freqs,noise); % zero frequencies, and zero delay

%% AR Solve
L_AR = opt_L(coeffs',alpha,beta);
A_AR = laplacian_to_adjacency(L_AR);
%% Dong Solve
L_Dong = full_optimization(y',alpha,beta,tolerance);
A_Dong = laplacian_to_adjacency(L_Dong);

%% Comparison
threshold_array = .01:.01:.3;
f = zeros(length(threshold_array),1);
for i = 1:length(threshold_array)
        t = threshold_array(i);
        f1(i) = calculate_F_score(A>0,A_AR>t);
        f2(i) = calculate_F_score(A>0,A_Dong>t);
end
figure
plot(f1)
figure
plot(f2)