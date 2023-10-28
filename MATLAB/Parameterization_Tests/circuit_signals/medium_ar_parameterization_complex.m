clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Signal Configuration
signal_params.N = 20;
signal_params.M = 1e3;
signal_params.active = 20;
signal_params.noise = .25;
signal_params.additional_diag = 2.5; %shunt resistors
signal_params.trace_normalization = signal_params.N;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;
signal_params.isComplex = 1;

%% GL_AR Configuration
AR_params.N = signal_params.N;
AR_params.P = 20;
AR_params.gamma = 500;%2.6;
AR_params.threshold = .35;%.085;
%% Parameterization
gammas = logspace(-8,8,20);
thresholds = 0:1e-5:.05;
noises = logspace(-4,0,10);
t_max = length(thresholds);
trials = 20;
F = zeros(length(gammas),length(thresholds),length(noises));
count = 0;
for j = 1:length(noises)
    for k = 1:length(gammas)
        toc
        tic
        count/(length(gammas)*length(noises))*100
        count = count+1;
        AR_params.gamma = gammas(k);
        signal_params.noise = noises(j);
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
        F(k,:,j) = mean(f,2);
    end
    max(F(:,:,j),[],"all")
end

%% Comparison
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
% surf(X_axis,Y_axis,F');
% s.EdgeColor = 'none';