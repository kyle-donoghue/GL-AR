clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Signal Configuration
signal_params.N = 20;
signal_params.M = 1e3;
signal_params.active = 20;
signal_params.noise = 1.8e-4;
signal_params.additional_diag = 8;
signal_params.trace_normalization = signal_params.N/60;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;


%% Dong Configuration
dong_params.N = signal_params.N;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);

%% Parameterization
betas = logspace(-8,-8,20);
thresholds = 0:1e-3:.15;
t_max = length(thresholds);
trials = 20;
F = zeros(length(betas),length(thresholds));
for k = 1:length(betas)
    toc
    tic
    k/length(betas)*100
    dong_params.beta = betas(k);
    f = zeros(trials,1);
    parfor i = 1:trials
        L_0 = graphs.create(signal_params,'er',0.2);
        V = signals.circuit_sine(signal_params,L_0);
        [L,~,~] = GL.dong(V,dong_params);
        for t = 1:t_max 
            L_tmp = GL.threshold(L,thresholds(t))
            [~,~,f(t,i),~,~] = graphs.performance(L_0,L_tmp);
        end
    end
    F(k,:) = mean(f,2);
end

%% Comparison
[X_axis,Y_axis] = meshgrid(log10(betas),thresholds);
surf(X_axis,Y_axis,F');