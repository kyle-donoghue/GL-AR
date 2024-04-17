clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Signal Configuration
signal_params.N = 20;
signal_params.M = 1e3;
signal_params.active = 20;
signal_params.noise = 1.8e-4;
signal_params.additional_diag = 8; %shunt resistors
signal_params.trace_normalization = signal_params.N/60;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;
signal_params.isComplex = 1;

%% GL_AR Configuration
AR_params.N = signal_params.N;
AR_params.P = 20;
AR_params.gamma = 2.6;
AR_params.threshold = .085;

%% Parameterization
gammas = logspace(0,16,20);
thresholds = 0:1e-3:.15;
t_max = length(thresholds);
trials = 20;

% noises = logspace(-8,0,10);
noises = 1.8e-4;
diags = logspace(-1,2,10);
traces = 10:10:100;

bigF = zeros(length(noises),length(diags),length(traces));

count = 0;
for noise_i = 1:length(noises)
    for diag_i = 1:length(diags)
        for trace_i = 1:length(traces)
            toc
            tic
            count/(length(noises)*length(diags)*length(traces))*100
            count = count+1;
            F = zeros(length(gammas),length(thresholds));
            for k = 1:length(gammas)

                AR_params.gamma = gammas(k);
                signal_params.noise = noises(noise_i);
                signal_params.additional_diag = diags(diag_i);
                signal_params.trace_normalization = signal_params.N/traces(trace_i);

                f = zeros(trials,1);
                parfor i = 1:trials
                    L_0 = graphs.create(signal_params,'er',0.2);
                    V = signals.circuit_sine(signal_params,L_0);
                    L = GL.AR(V,AR_params);
                    for t = 1:t_max 
                        L_tmp = GL.threshold(L,thresholds(t))
                        [~,~,f(t,i),~,~] = graphs.performance(L_0,L_tmp);
                    end
                end
                F(k,:) = mean(f,2);
            end
            bigF(noise_i,diag_i,trace_i) = max(F,[],"all");
        end
    end
end

%% Comparison
max(bigF,[],"all")