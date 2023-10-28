clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
addpath('../progressBar/');
%% INITIALIZATION
n = 20;
noise = .25;
N = 1e3;
delay = 50;
interval = 5;
Fs = 250;
P = 30;
ER_prob = .2;

alpha = 1;
beta = .15;

threshold = .065;
%% STATISTICAL STORAGE
trials = 5e2;
f_scores = zeros(trials,1);
ppm = ParforProgressbar(trials);
%% SIMULATION AND OPTIMIZATION
parfor k = 1:trials
    ppm.increment();
    
    freqs = generate_frequencies2(n,5,100);
    [coeffs, y, x, A] = sineER_random(n,P,ER_prob,delay,interval,N,Fs,freqs,noise);
    L_interval = zeros(n,n,interval);
    A_interval = zeros(n,n,interval);
    for i = 1:interval
        [L_interval(:,:,i)] = opt_L(coeffs(:,:,i)',alpha,beta);
        A_interval(:,:,i) = laplacian_to_adjacency(L_interval(:,:,i));
    end
    % post-processing
    
    A_new = mean(A_interval>threshold,3);
    f_scores(k) = calculate_F_score(A>0,A_new);
end
delete(ppm);

%% FIGURES
figure
histogram(f_scores);