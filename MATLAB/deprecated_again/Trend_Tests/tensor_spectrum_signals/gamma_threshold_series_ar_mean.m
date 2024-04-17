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
for i = 1:8
    s_p(i) = signal_params;
end
s_p(2).minSeparation = 50;
s_p(2).maxSeparation = 250;
s_p(3).minSeparation = 5;
s_p(3).maxSeparation = 50;
s_p(4).min_delay = 0;
s_p(4).max_delay = 20;
s_p(5).min_delay = 50;
s_p(5).max_delay = 100;
s_p(6).interval_length = 2048;

for i = 1:length(s_p)
    s_p(i).M = s_p(i).intervals*s_p(i).interval_length;
    a_p(i) = s_p(i);
end
for i = 1:8
    a_p(i).P = 40;
    a_p(i).gamma = 0;
    a_p(i).threshold = 0;
end
a_p(7).P = 80;
a_p(8).P = 20;

%%
AR_params = signal_params;
AR_params.P = 40;
AR_params.gamma = 0;
AR_params.threshold = 0;
%% Parameter Series
gammas = logspace(.5,2,40);
thresholds = 0:.005:.3;
t_max = length(thresholds);
trials = 20;

indice_mat = zeros(length(s_p),3);
F_mat = zeros(length(gammas),length(thresholds),length(s_p));

count = 0;
for s = 1:length(s_p)
    F = zeros(length(gammas),length(thresholds));
    signal_params = s_p(s);
    AR_params = a_p(s);
    for k = 1:length(gammas)
        toc
        tic
        count/(length(s_p)*length(gammas))*100
        count = count+1;
        AR_params.gamma = gammas(k);
        f = zeros(trials,1);
        parfor i = 1:trials
            % [L_0,A] = graphs.create(signal_params,'er',.2);
            % [L_0,A] = graphs.create(signal_params,'pa',1);
            [L_0,A] = graphs.create(signal_params,'gaussian',0.75,0.5);
            A_0 = graphs.to_directed(A);
            G = graphs.createGraphTensor(signal_params,A_0);
            x = signals.create_raw_sine(signal_params);
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
    F_mat(:,:,s) = F;
    [maximum, ind] = max(F,[],"all");
    [row, col] = ind2sub(size(F),ind);
    indice_mat(s,1) = gammas(row);
    indice_mat(s,2) = thresholds(col);
    indice_mat(s,3) = maximum;    
end

