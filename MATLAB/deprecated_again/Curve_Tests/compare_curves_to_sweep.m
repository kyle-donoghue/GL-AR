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
AR_params.threshold = 0;
%% Parameterization
sparsity = .1;
gammas = logspace(-2,4,20);
thresholds = 0:.005:.5;
t_max = length(thresholds);
trials = 10;
F = zeros(length(gammas),length(thresholds));
for k = 1:length(gammas)
    toc
    tic
    k/length(gammas)*100
    AR_params.gamma = gammas(k);
    f = zeros(trials,1);
    parfor i = 1:trials
        [L_0,A] = graphs.create(signal_params,'er',sparsity);
        % [L_0,A] = graphs.create(signal_params,'pa',1);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        x = signals.createFullRectPulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        L = GL.AR_mean(y_noisy,AR_params);

        for t = 1:t_max 
            L_tmp = GL.threshold(L,thresholds(t));
            [~,~,f(t,i),~,~] = graphs.performance(L_0,L_tmp);
        end
    end
    F(k,:) = mean(f,2);
end

%% Comparison
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
s = surf(X_axis,Y_axis,F');hold on;
s.EdgeColor = "none";
%%
log_b = -7.177;
log_c = 0.2361;
log_d = 1.806;
log_offset = -.12;
exp_a = 0.03537;
exp_b = 1.941;
exp_c = .95*.05;

AR_params.gamma = find_gamma(sparsity,log_b,log_c,log_d,log_offset)
AR_params.threshold = find_threshold(AR_params.gamma,exp_a,exp_b,exp_c)

%%
scatter3(log10(AR_params.gamma),AR_params.threshold,max(F,[],"all"),'filled','MarkerFaceColor','red')
%%
function t = find_threshold(x,a,b,c)
    t = a*exp(b*(log10(x)-1))+c;
end
function g = find_gamma(x,b,c,d,offset)
    tmp = d + (-d)/(1 + (x/c)^b)+offset; % the output is in log
    g = 10^tmp;
end