clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Signal Configuration
signal_params.N = 20;
signal_params.M = 1e3;
signal_params.active = 20;
signal_params.SNR = 2;
signal_params.additional_diag = .15; %shunt resistors
signal_params.trace_normalization = signal_params.N;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;
signal_params.isComplex = 0;

%% GL_AR Configuration
AR_params.N = signal_params.N;
AR_params.P = 20;
AR_params.gamma = 500;%2.6;
AR_params.threshold = .35;%.085;
%% Parameterization
gammas = logspace(-4,4,20);
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
        L_0 = graphs.create(signal_params,'er',0.2);
        V = signals.circuit_sine(signal_params,L_0);
        L = GL.AR(V,AR_params);
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