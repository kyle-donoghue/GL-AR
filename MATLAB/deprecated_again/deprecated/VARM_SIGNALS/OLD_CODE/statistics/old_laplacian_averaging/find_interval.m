clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
addpath('../progressBar/');
D = parallel.pool.DataQueue;
afterEach(D, @nUpdateWaitbar);
clear nUpdateWaitbar;
%% INITIALIZATION
n = 20;
noise = .25;
N = 1e3;
delay = 50;
Fs = 250;
% freqs = generate_frequencies(n,1,100,.1);
ER_prob = .2;
P = 30;



%% STATISTICAL STORAGE
trials = 1e2;

alpha = 1;
beta = .15;

intervals = 1:1:30;

f_scores_final = zeros(length(intervals),2);
%% SIMULATION AND OPTIMIZATION
for b = 1:length(intervals)
    interval = intervals(b);
    b
    nUpdateWaitbar(length(intervals));
    f_scores_tmp = zeros(trials,1);
    parfor k = 1:trials
        freqs = generate_frequencies2(n,5,100);
        [coeffs, y, x, A] = sineER(n,P,ER_prob,delay,interval,N,Fs,freqs,noise);
        L_interval = zeros(n,n,interval);
        for i = 1:interval
            [L_interval(:,:,i)] = opt_L(coeffs(:,:,i)',alpha,beta);
        end
        L = mean(L_interval,3);
    
        % post-processing
        A_new = laplacian_to_adjacency(L);
        threshold_array = .01:.01:.3;
        f = zeros(length(threshold_array),1);
        for i = 1:length(threshold_array)
            t = threshold_array(i);
            f(i) = calculate_F_score(A>0,A_new>t);
        end
        f(isnan(f)) = 0;
        f_scores_tmp(k) = max(f);
    end
    f_scores_final(b,:) = [mean(f_scores_tmp) std(f_scores_tmp)];
end
%% DATA STORAGE
save('interval_sim.mat')