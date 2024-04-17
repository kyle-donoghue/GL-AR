clc;clear;close all hidden;
addpath(genpath('GL_classes/'));
%%
trials = 32;
er_param = .2;
%% ER=.2
%N_list = 100:-4:44;
N_list = 68:-4:12;
N_list = 10;

F_list = zeros(length(N_list),1);
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,'er');
    signal_params.SNR = .02;
    [F_list(i),F_surfs(:,:,i)] = run_maximum_sims(signal_params,er_param,trials);
    i/length(N_list)*100
    toc
    tic
end
%%
save("large_maximum_data_results_count_er_signaltest_SNR02.mat")
% quit
%% 
function [max_F,sparse_F] = run_maximum_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);

    % AR_params.P = 100;
    
    thresholds = 0:.0025:2;
    gammas = logspace(-2,4,60);
    t_max = length(thresholds);

    sparse_F = zeros(length(gammas),length(thresholds),length(graph_params));
    sparsities = zeros(length(graph_params),1);
    count = 1;
    for s = 1:length(graph_params)
        p = graph_params(s);
        sparsities(s) = graphs.get_sparsity(signal_params.graph_type,p);
        F = zeros(length(gammas),length(thresholds));
        for k = 1:length(gammas)
            fprintf("%d %s",signal_params.N,signal_params.graph_type)
            tim = toc
            tic
            percent = count/(length(gammas)*length(graph_params))*100
            left = tim*((length(gammas)*length(graph_params))-count)/60
            count = count+1;
    
            
            AR_params.gamma = gammas(k);
        
            f = zeros(t_max,trials);
            parfor i = 1:trials
                [L_0,~,A_d] = graphs.create(signal_params,p);
                G = graphs.createGraphTensor(signal_params,A_d);
                y_noisy = signals.generateFilteredRectPulse(signal_params,G);
            
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
    max_F = max(sparse_F,[],"all");
end
