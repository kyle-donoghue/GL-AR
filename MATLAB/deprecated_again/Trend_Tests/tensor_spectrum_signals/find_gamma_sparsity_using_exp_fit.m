clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Fit Data (y = a*exp(b*(log10(x)-1))+c)
a = 0.03537;
b = 1.941;
c = .95*.05;
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
AR_params.threshold = 0;
%%
gammas = logspace(-1.5,3,40);
sparsities = .05:.025:.95;
trials = 20;
max_gamma = zeros(length(sparsities),1);
gamma_curves = zeros(length(sparsities),length(gammas));
count = 1;
for s = 1:length(sparsities)
    F = zeros(length(gammas),1);
    for k = 1:length(gammas)
        tim = toc
        tic
        percent = count/(length(gammas)*length(sparsities))*100
        left = tim*((length(gammas)*length(sparsities))-count)/60
        count = count+1;
        AR_params.gamma = gammas(k);
        f = zeros(trials,1);
        parfor i = 1:trials
            [L_0,A] = graphs.create(signal_params,'er',sparsities(s));
            % [L_0,A] = graphs.create(signal_params,'pa',1);
            A_0 = graphs.to_directed(A);
            G = graphs.createGraphTensor(signal_params,A_0);
            x = signals.createFullRectPulse(signal_params);
            X = signals.createTensorSpectrum(signal_params,x);
            Y = signals.filterTensorSpectrum(signal_params,X,G);
            y = signals.inverseTensorSpectrum(signal_params,Y);
            y_noisy = signals.add_noise(signal_params,y);
            
            L = GL.AR_mean(y_noisy,AR_params);
            
            L_tmp = GL.threshold(L,find_threshold(AR_params.gamma,a,b,c));
            [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
        end
        F(k) = mean(f);
    end
    gamma_curves(s,:) = F;
    max_gamma(s) = gammas(F == max(F));
end
save("gamma_sparsity_relationship_v3.mat")
function t = find_threshold(x,a,b,c)
    t = a*exp(b*(log10(x)-.75))+c;
end