clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
signal_params.N = 20;
signal_params.raw = 1;

signal_params.SNR = 2e0;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;

signal_params.intervals = signal_params.N;
signal_params.interval_length = 512;
signal_params.M = signal_params.intervals*signal_params.interval_length;

signal_params.order = 1;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .1;
signal_params.pole_mean = .5;

signal_params.min_delay = 10;
signal_params.max_delay = 50;


signal_params.minSeparation = 100;
signal_params.maxSeparation = 500;

%%
AR_params = signal_params;
AR_params.P = 40;
AR_params.gamma = 0;
AR_params.threshold = .1;

%% Parameter Series
gammas = logspace(-2,4,40);
thresholds = .005:.005:.3;
t_max = length(thresholds);

[L_0,A] = graphs.create(signal_params,'er',.2);
% [L_0,A] = graphs.create(signal_params,'pa',1);
% [L_0,A] = graphs.create(signal_params,'gaussian',0.75,0.5);

A_0 = graphs.to_directed(A);
G = graphs.createGraphTensor(signal_params,A_0);
x = signals.create_raw_sine(signal_params);
X = signals.createTensorSpectrum(signal_params,x);
Y = signals.filterTensorSpectrum(signal_params,X,G);
y = signals.inverseTensorSpectrum(signal_params,Y);
y_noisy = signals.add_noise(signal_params,y);

trials = 500;

weight_size = 1/2*signal_params.N*(signal_params.N-1);

f = zeros(length(gammas),t_max);
F = f;
weights = zeros(length(gammas),weight_size);
big_weights = zeros(length(gammas),weight_size*trials);

counter = 1;
for t = 1:trials
    time = toc
    tic
    percent = counter/(trials)*100
    time_left = ((trials)-counter)*time/60
    counter = counter+1;
    parfor i = 1:length(gammas)
        ar_p = AR_params;
        ar_p.gamma = gammas(i);
    
        L = GL.AR_mean(y_noisy,ar_p);
    
        for thresh = 1:t_max 
            L_tmp = GL.threshold(L,thresholds(thresh));
            [~,~,f(i,thresh),~,~] = graphs.performance(L_0,L_tmp);
        end
        A_tmp = graphs.to_adjacency(L);
        weights(i,:) = squareform((A_tmp-diag(A_tmp)).');
    end
    F = F+f;
    big_weights(:,(t-1)*weight_size+1:(t)*weight_size) = weights;
end
F = F/trials;

[~,~,~,~,num_edge] = graphs.performance(L_0,L_0);
sparsity = num_edge/(signal_params.N*(signal_params.N-1)/2)

%%
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
s = surf(X_axis,Y_axis,f');
s.EdgeColor = "none";

%%
n_hist = 25;
figure
subplot(1,3,1)
histogram(big_weights(20,:),n_hist);
xline(thresholds(12));
title(sprintf("Gamma=%.1e, F=%.1d",gammas(20),f(20,12)*100));
subplot(1,3,2)
histogram(big_weights(10,:),n_hist);
xline(thresholds(12));
title(sprintf("Gamma=%.1e, F=%.1d",gammas(15),f(10,12)*100));
subplot(1,3,3)
histogram(big_weights(3,:),n_hist);
xline(thresholds(12));
title(sprintf("Gamma=%.1e, F=%.1d",gammas(3),f(3,12)*100));

%%
figure;
for i = length(gammas):-1:1
    histogram(big_weights(i,:),n_hist);
    th = .1;
    xline(th);
    xlim([0 1]);
    title(sprintf("Gamma=%d, F=%d",gammas(i),f(i,find(thresholds==th))));
    pause(.5)
end
