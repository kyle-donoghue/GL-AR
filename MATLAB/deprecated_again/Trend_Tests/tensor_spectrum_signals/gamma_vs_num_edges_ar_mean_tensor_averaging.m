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
AR_params = signal_params;
AR_params.P = 40;
AR_params.gamma = 0;
AR_params.threshold = .1;
%% Parameter Series
gammas = logspace(0,2,40);
thresholds = logspace(-3,-.5,100);
t_max = length(thresholds);

graph_params = .05:.05:.8;

trials = 50;

scatters = zeros(trials*length(graph_params),4);

counter = 1;
for k = 1:length(graph_params)
    for t = 1:trials
        time = toc
        tic
        percent = counter/(trials*length(graph_params))*100
        time_left = ((trials*length(graph_params))-counter)*time/60

        [L_0,A] = graphs.create(signal_params,'er',graph_params(k));
        % [L_0,A] = graphs.create(signal_params,'pa',1);
        % [L_0,A] = graphs.create(signal_params,'gaussian',0.75,0.5);
        A_0 = graphs.to_directed(A);
        G = graphs.createGraphTensor(signal_params,A_0);
        x = signals.create_raw_sine(signal_params);
        X = signals.createTensorSpectrum(signal_params,x);
        Y = signals.filterTensorSpectrum(signal_params,X,G);
        y = signals.inverseTensorSpectrum(signal_params,Y);
        y_noisy = signals.add_noise(signal_params,y);
    
        f = zeros(length(gammas),t_max);
        parfor i = 1:length(gammas)
            ar_p = AR_params;
            ar_p.gamma = gammas(i);
            
            L = GL.AR_mean(y_noisy,ar_p);

            for thresh = 1:t_max 
                L_tmp = GL.threshold(L,thresholds(thresh));
                [~,~,f(i,thresh),~,~] = graphs.performance(L_0,L_tmp);
            end
        end
        [max_f, ind] = max(f,[],"all");
        [row, col] = ind2sub(size(f),ind);
        max_gamma = gammas(row);
        max_threshold = thresholds(col);
        
        [~,~,~,~,num_edge] = graphs.performance(L_0,L_0);
        sparsity = num_edge/(signal_params.N*(signal_params.N-1)/2);
        
        scatters(counter,:) = [sparsity max_gamma max_threshold max_f];

        counter=counter+1;
    end
end
%%
figure;
subplot(1,2,1);
scatter3(scatters(:,1),scatters(:,2),scatters(:,4),'filled');
xlabel('Sparsity')
ylabel('Gamma')
subplot(1,2,2);
scatter3(scatters(:,1),scatters(:,3),scatters(:,4),'filled');
xlabel('Sparsity')
ylabel('Threshold')
