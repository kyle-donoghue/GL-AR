clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
addpath('../progressBar/');
%% INITIALIZATION
n = 20;
noise = .25;
N = 1e3;
delay = 1;
interval = 5;
Fs = 250;
% freqs = generate_frequencies(n,1,100,.1);
freqs = generate_frequencies2(n,5,100);
P = 30;
ER_prob = .2;

alpha = 1;
beta = .15;

A = adjacencyER(n,ER_prob);

threshold = .065;
%% STATISTICAL STORAGE
trials = 1e2;
f_scores = zeros(trials,1);
ppm = ParforProgressbar(trials);
%% SIMULATION AND OPTIMIZATION
parfor k = 1:trials
    ppm.increment();
    
    x = createSine(n,N,Fs,freqs,noise);
    interval_length = round(N/interval);
    y = zeros(interval_length,n,interval);
    coeffs = zeros(P+1,n,interval);
    if delay == 0
        D_AR = eye(n)+A;
    else
        mdl = createVARM(n, delay, A);
        D_AR = cat(3,mdl.AR{:});
        D_AR(:,:,2:end+1) = D_AR(:,:,1:end);
        D_AR(:,:,1) = eye(n);
    end
    for i = 1:interval
        y(:,:,i) = customFilter(D_AR,x(((i-1)*interval_length+1):(i*interval_length),:),n);
        coeffs(:,:,i) = approxAR(y(:,:,i),P);
    end
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
mean(f_scores)