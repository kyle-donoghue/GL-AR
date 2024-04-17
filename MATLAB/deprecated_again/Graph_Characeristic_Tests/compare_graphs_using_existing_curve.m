clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
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

trials = 8;
%%
AR_params = signal_params;
AR_params.P = 40;
AR_params.gamma = 0;
AR_params.threshold = 0;
%% ER Simulation
er_param = .05:.05:.95;
F_ER = zeros(length(er_param),1);
sparsities_ER = zeros(length(er_param),1);
for s = 1:length(er_param)
    toc
    tic
    s/length(er_param)*100
    p = er_param(s);

    sparsity = graphs.get_sparsity('er',p);
    edge_spread = graphs.get_edge_spread('er',p,signal_params.N);
    
    AR_params.threshold = GL.get_threshold(edge_spread,signal_params.N);
    AR_params.gamma = GL.get_gamma(sparsity,AR_params.threshold,signal_params.N);
    
    sparsities_ER(s) = sparsity;
    
    f = zeros(trials,1);
    parfor i = 1:trials
        [L_0,A] = graphs.create(signal_params,'er',p);
        % [L_0,A] = graphs.create(signal_params,'pa',M(s));
        % [L_0,A] = graphs.create(signal_params,'gaussian',0.75,M(s));
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        x = signals.createFullRectPulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        L = GL.AR_mean(y_noisy,AR_params);

        L_tmp = GL.threshold(L,AR_params.threshold);
        [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
    end
    F_ER(s) = mean(f);
end
%% ER Comparison Against Maximum
load('maximums_F.mat')
figure;hold on;
plot(max_sparsities_ER,maximums_ER,'--');
plot(sparsities_ER,F_ER*100)
title("ER")

%% BA Simulation
ba_param = 1:10;
F_BA = zeros(length(ba_param),1);
sparsities_BA = zeros(length(ba_param),1);
for s = 1:length(ba_param)
    toc
    tic
    s/length(ba_param)*100
    p = ba_param(s);

    sparsity = graphs.get_sparsity('ba',p);
    edge_spread = graphs.get_edge_spread('ba',p,signal_params.N);
    
    AR_params.threshold = GL.get_threshold(edge_spread,signal_params.N);
    AR_params.gamma = GL.get_gamma(sparsity,AR_params.threshold,signal_params.N);
    
    sparsities_BA(s) = sparsity;
    
    f = zeros(trials,1);
    parfor i = 1:trials
        [L_0,A] = graphs.create(signal_params,'pa',p);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        x = signals.createFullRectPulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        L = GL.AR_mean(y_noisy,AR_params);

        L_tmp = GL.threshold(L,AR_params.threshold);
        [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
    end
    F_BA(s) = mean(f);
end
%% BA Comparison Against Maximum
load('maximums_F.mat')
figure;hold on;
plot(max_sparsities_BA,maximums_BA,'--');
plot(sparsities_BA,F_BA*100)
title("BA")

%% RND Simulation
rnd_param = .1:.05:1.2;
F_RND = zeros(length(rnd_param),1);
sparsities_RND = zeros(length(rnd_param),1);
for s = 1:length(rnd_param)
    toc
    tic
    s/length(rnd_param)*100
    p = rnd_param(s);

    sparsity = graphs.get_sparsity('gaussian',p);
    edge_spread = graphs.get_edge_spread('gaussian',p,signal_params.N);
    
    AR_params.threshold = GL.get_threshold(edge_spread,signal_params.N);
    AR_params.gamma = GL.get_gamma(sparsity,AR_params.threshold,signal_params.N);
    
    sparsities_RND(s) = sparsity;
    
    f = zeros(trials,1);
    parfor i = 1:trials
        [L_0,A] = graphs.create(signal_params,'gaussian',0.75,p);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        x = signals.createFullRectPulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        L = GL.AR_mean(y_noisy,AR_params);

        L_tmp = GL.threshold(L,AR_params.threshold);
        [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
    end
    F_RND(s) = mean(f);
end
%% RND Comparison Against Maximum
load('maximums_F.mat')
figure;hold on;
plot(max_sparsities_RND,maximums_RND,'--');
plot(sparsities_RND,F_RND*100)
title("RND")