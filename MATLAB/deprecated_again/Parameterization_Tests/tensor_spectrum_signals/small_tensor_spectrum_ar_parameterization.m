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

%%
AR_params.N = signal_params.N;
AR_params.P = 10;
AR_params.gamma = 0;
AR_params.threshold = 0;
%% Parameterization
gammas = logspace(.5,2,40);
thresholds = logspace(-4,-1,25);
t_max = length(thresholds);
trials = 20;
F = zeros(length(gammas),length(thresholds));
for k = 1:length(gammas)
    toc
    tic
    k/length(gammas)*100
    AR_params.gamma = gammas(k);
    f = zeros(trials,1);
    parfor i = 1:trials
        [L_0,A] = graphs.create(signal_params,'er',.2);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        % x = signals.create_raw_sine(signal_params);
        x = signals.createFullRectPulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        L = GL.AR(y_noisy,AR_params);
        for t = 1:t_max 
            L_tmp = GL.threshold(L,thresholds(t));
            [~,~,f(t,i),~,~] = graphs.performance(L_0,L_tmp);
        end
    end
    F(k,:) = mean(f,2);
end

%% Comparison
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
s = surf(X_axis,Y_axis,F');
s.EdgeColor = "none";

%%
[maximum, ind] = max(F,[],"all");
[row, col] = ind2sub(size(F),ind);
gamma = gammas(row)
thresh = thresholds(col)
maximum