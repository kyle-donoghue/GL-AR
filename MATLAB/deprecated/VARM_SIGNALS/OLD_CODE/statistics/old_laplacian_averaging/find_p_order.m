clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
addpath('../progressBar/');
D = parallel.pool.DataQueue;
afterEach(D, @nUpdateWaitbar);
%% INITIALIZATION
n = 20;
noise = .25;
N = 1e3;
delay = 50;
interval = 5;
Fs = 250;
% freqs = generate_frequencies(n,1,100,.1);
freqs = generate_frequencies2(n,5,100);
ER_prob = .2;



A = adjacencyER(n,ER_prob);
mdl = createVARM(n, delay, A);

%% STATISTICAL STORAGE
trials = 1e2;

alpha = 1;
beta = .15;

orders = 10:2:100;

f_scores_final = zeros(length(orders),2);
%% SIMULATION AND OPTIMIZATION
for b = 1:length(orders)
    P = orders(b);
    b
    nUpdateWaitbar(length(orders));
    f_scores_tmp = zeros(trials,1);
    parfor k = 1:trials
        
        x = createSine(n,N,Fs,freqs,noise);
        interval_length = round(N/interval);
        y = zeros(interval_length,n,interval);
        coeffs = zeros(P+1,n,interval);
        D_AR = cat(3,mdl.AR{:});
        D_AR(:,:,2:end+1) = D_AR(:,:,1:end);
        D_AR(:,:,1) = eye(n);
        for i = 1:interval
            y(:,:,i) = customFilter(D_AR,x(((i-1)*interval_length+1):(i*interval_length),:),n);
            coeffs(:,:,i) = approxAR(y(:,:,i),P);
        end
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
save('p_order_sim.mat')