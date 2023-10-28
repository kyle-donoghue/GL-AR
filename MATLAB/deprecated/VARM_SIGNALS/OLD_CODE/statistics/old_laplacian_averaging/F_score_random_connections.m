clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
addpath('../progressBar/');
%% INITIALIZATION
n = 20;
noise = .1;
N = 1e3;
delay = 50;
interval = 5;
Fs = 250;
% freqs = generate_frequencies(n,1,100,.1);
freqs = generate_frequencies2(n,5,100);
P = 30;
ER_prob = .2;

alpha = 1;
beta = .15;

%% STATISTICAL STORAGE
trials = 5e2;
f_scores = zeros(trials,1);
ppm = ParforProgressbar(trials);
%% SIMULATION AND OPTIMIZATION
parfor k = 1:trials
    ppm.increment();
    [coeffs, y, x, A] = sineER(n,P,ER_prob,delay,interval,N,Fs,freqs,noise);
    L_interval = zeros(n,n,interval);
    for i = 1:interval
        [L_interval(:,:,i)] = opt_L(coeffs(:,:,i)',alpha,beta);
    end
    L = mean(L_interval,3);

    % post-processing
    A = A>0;
    A_new = laplacian_to_adjacency(L);
    threshold_array = .01:.01:.3;
    f = zeros(length(threshold_array),1);
    for i = 1:length(threshold_array)
        t = threshold_array(i);
        f(i) = calculate_F_score(A>0,A_new>t);
    end
    f(isnan(f)) = 0;
    f_scores(k) = max(f);
end
delete(ppm);

%% FIGURES
histogram(f_scores);