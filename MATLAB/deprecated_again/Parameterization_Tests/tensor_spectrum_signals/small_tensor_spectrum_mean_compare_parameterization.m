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

signal_params.order = 5;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .1;
signal_params.pole_mean = .5;

signal_params.min_delay = 10;
signal_params.max_delay = 50;


signal_params.minSeparation = 100;
signal_params.maxSeparation = 500;
% signal_params.minSeparation = 50;
% signal_params.maxSeparation = 200;

%%
AR_params = signal_params;
AR_params.P = 40;
AR_params.gamma = 0;
AR_params.threshold = 0;
%%
dong_params.N = signal_params.N;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
%% Parameterization
gammas = logspace(-3,3,20);
thresholds = 0:.001:.2;
t_max = length(thresholds);
trials = 20;
F_AR = zeros(length(gammas),length(thresholds));
F_Dong = zeros(length(gammas),length(thresholds));
for k = 1:length(gammas)
    toc
    tic
    k/length(gammas)*100
    AR_params.gamma = gammas(k);
    dong_params.beta = gammas(k);
    f_AR = zeros(trials,1);
    f_Dong = zeros(trials,1);
    parfor i = 1:trials
        % [L_0,A] = graphs.create(signal_params,'er',.2);
        % [L_0,A] = graphs.create(signal_params,'pa',1);
        [L_0,A] = graphs.create(signal_params,'gaussian',0.75,0.5);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);

        x = signals.createFullRectPulse(signal_params);
        % x = signals.create_raw_sine(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);

        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        L_AR = GL.AR_mean(y_noisy,AR_params);
        L_Dong = GL.dong(y_noisy,dong_params);

        for t = 1:t_max 
            L_AR_tmp = GL.threshold(L_AR,thresholds(t));
            L_Dong_tmp = GL.threshold(L_Dong,thresholds(t));
            [~,~,f_AR(t,i),~,~] = graphs.performance(L_0,L_AR_tmp);
            [~,~,f_Dong(t,i),~,~] = graphs.performance(L_0,L_Dong_tmp);
        end
    end
    F_AR(k,:) = mean(f_AR,2);
    F_Dong(k,:) = mean(f_Dong,2);
end

%% Comparison
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
figure;
subplot(1,2,1);
s = surf(X_axis,Y_axis,F_AR');
s.EdgeColor = "none";
subplot(1,2,2);
s = surf(X_axis,Y_axis,F_Dong');
s.EdgeColor = "none";