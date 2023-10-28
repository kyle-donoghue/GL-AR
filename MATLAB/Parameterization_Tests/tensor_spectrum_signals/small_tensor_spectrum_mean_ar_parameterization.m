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
signal_params.pole_mean = .7;

signal_params.minSeparation = 20;
signal_params.maxSeparation = 200;

%%
AR_params.N = signal_params.N;
AR_params.P = 30;
AR_params.gamma = 0;
AR_params.threshold = 0;
%% Parameterization
gammas = logspace(-3,3,20);
thresholds = 0:.001:.2;
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
        x = signals.createFullSinePulse(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
        
        B = zeros(signal_params.N,AR_params.P,signal_params.intervals);
        for j = 1:signal_params.intervals
            B(:,:,j) = signals.approxAR(y_noisy(:,(1+(j-1)*signal_params.interval_length):(j*signal_params.interval_length)),AR_params.P);
        end
        B = mean(B,3);

        B = B-mean(B,2);

        Q = create_Q_matrix(AR_params.N,1);
        c = create_c_vec(AR_params.N,B)*AR_params.gamma;
        A = create_constraint_matrix(AR_params.N);
        b = create_constraint_vec(AR_params.N,AR_params.N);

        lowerbound = -Inf(size(b));
        lowerbound(1:AR_params.N+1) = b(1:AR_params.N+1);
        upperbound = b;
        m = osqp;
        m.setup(Q,c,A,lowerbound,upperbound,'verbose',false);
        results = m.solve();
        phi = results.x;
        
        l = create_dup_matrix(AR_params.N)*phi;
        L = convert_to_matrix(l);

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

