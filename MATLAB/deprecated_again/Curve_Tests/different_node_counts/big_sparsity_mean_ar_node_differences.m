clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
signal_params.N = 10;
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
% sparsities = .1:.05:.95;
M = 1:12;
% M = .1:.1:1.2;
sparsities = zeros(length(M),1);

gammas = logspace(-2,4,40);
thresholds = 0:.0025:2;
t_max = length(thresholds);
trials = 8;

sparse_F = zeros(length(gammas),length(thresholds),length(sparsities));
count = 1;
for s = 1:length(M)
    F = zeros(length(gammas),length(thresholds));
    sparsities(s) = graphs.get_sparsity('pa',M(s));
    for k = 1:length(gammas)
        tim = toc
        tic
        percent = count/(length(gammas)*length(M))*100
        left = tim*((length(gammas)*length(M))-count)/60
        count = count+1;
        AR_params.gamma = gammas(k);
        f = zeros(t_max,trials);
        parfor i = 1:trials
            % [L_0,A] = graphs.create(signal_params,'er',sparsities(s));
            [L_0,A] = graphs.create(signal_params,'pa',M(s));
            % [L_0,A] = graphs.create(signal_params,'gaussian',0.75,M(s));
            while (~any(L_0,'all') || all(L_0,'all'))
                [L_0,A] = graphs.create(signal_params,'pa',M(s));
            end
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
    sparse_F(:,:,s) = F;
end
save("gamma_threshold_sweep_with_sparsity_big_n10_ba.mat")
