clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
load("maximums_sweep.mat")
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
sparsities = .05:.025:.95;
trials = 128;

calculated = zeros(length(sparsities),1);
count = 1;
for s = 1:length(sparsities)
    tim = toc
    tic
    percent = count/(length(sparsities))*100
    left = tim*(length(sparsities)-count)/60
    count = count+1;

    
    [AR_params.gamma, AR_params.threshold] = find_gamma_threshold(sparsities(s),signal_params.N);



    f = zeros(trials,1);
    parfor i = 1:trials
        [L_0,A] = graphs.create(signal_params,'er',sparsities(s));
        % [L_0,A] = graphs.create(signal_params,'pa',1);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        x = signals.createFullRectPulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        L = GL.AR_mean(y_noisy,AR_params);

        L_tmp = GL.threshold(L,AR_params.threshold);
        [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
    end
    calculated(s) = mean(f)*100;
end
save("calculated_curves.mat",'calculated')
%%
figure;hold on;
plot(sparsities,maximums,'--');
plot(sparsities,calculated)
ylim([0 100])
xlabel('Sparsity')
ylabel('F-Score')

%%
function [g,t] = find_gamma_threshold(s,N)
    a_a = 0.1289;
    a_b = -5.855;
    a_offset = 0.0060;
    b_a1 =  0.4011;
    b_c1 = 0.2049;
    b_offset = 1.5;
    cliff_curve_offset = 1/N;
    cliff_curve_delay =1;
    t_a = 0.4285;
    t_b = -0.06768;
    t_offset = .75/N;
    t_initial = 9.5;

    % convert sparsity to number of edges
    num_edges = s*(signal_params.N*(signal_params.N-1)/2);

    % calculate threshold using num_edges
    t = t_a*exp(t_b*(num_edges-t_initial))+t_offset;

    % calculate parameters of current cliff curve
    a = a_a*exp(a_b*s)+a_offset;
    b = b_a1*exp(-((s)/b_c1)^2)+b_offset;

    % calculate gamma based on current cliff curve
    if t > cliff_curve_offset/2
        g = max((log((t-.75*cliff_curve_offset)/a)/b+cliff_curve_delay),0);
    else
        g = 0;
    end
    g = 10^g;
end