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

%%
AR_params = signal_params;
AR_params.P = 40;
AR_params.gamma = 0;
AR_params.threshold = 0;
%%
log_b = -3.996;
log_c = 0.2361;
log_d = 1.806;
log_offset = -.4;
exp_a = 0.03537;
exp_b = 1.941;
exp_c = .95*.05;

%% Parameter Series
trials = 20;

sparsity = .2;
AR_params.gamma = find_gamma(sparsity,log_b,log_c,log_d,log_offset);
AR_params.threshold = find_threshold(AR_params.gamma,exp_a,exp_b,exp_c)
f = zeros(trials,1);
parfor i = 1:trials
    [L_0,A] = graphs.create(signal_params,'er',sparsity);
    % [L_0,A] = graphs.create(signal_params,'pa',1);
    % [L_0,A] = graphs.create(signal_params,'gaussian',0.75,0.5);
    A_0 = graphs.to_directed(A);
    G = graphs.createGraphTensor(signal_params,A_0);
    x = signals.create_raw_sine(signal_params);
    X = signals.createTensorSpectrum(signal_params,x);
    Y = signals.filterTensorSpectrum(signal_params,X,G);
    y = signals.inverseTensorSpectrum(signal_params,Y);
    y_noisy = signals.add_noise(signal_params,y);
    
    L = GL.AR_mean(y_noisy,AR_params);
    L_tmp = GL.threshold(L,AR_params.threshold);
    [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
end
F = mean(f)

function t = find_threshold(x,a,b,c)
    t = a*exp(b*(log10(x)-1))+c;
end
function g = find_gamma(x,b,c,d,offset)
    tmp = d + (-d)/(1 + (x/c)^b)+offset; % the output is in log
    g = 10^tmp;
end