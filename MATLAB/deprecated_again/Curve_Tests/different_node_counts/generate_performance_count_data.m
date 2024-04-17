clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 128;
rnd_param = .4;
%% ER=.2
N_list = 40:-2:6;
F_list = zeros(length(N_list),1);
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,'gaussian');
    [F_list(i),~] = run_performance_sims(signal_params,rnd_param,trials);
    i/length(N_list)*100
    toc
    tic
end
%%
save("large_performance_data_results_diff_n.mat")
%% 
function [F,sparsities] = run_performance_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);
    F = zeros(length(graph_params),1);
    sparsities = zeros(length(graph_params),1);
    for s = 1:length(graph_params)
    
        p = graph_params(s);
        sparsities(s) = graphs.get_sparsity(signal_params.graph_type,p);
        edge_spread = graphs.get_edge_spread2(signal_params.graph_type,p,signal_params.N);
        vertice_std = graphs.get_vertice_std(signal_params.graph_type,p,signal_params.N);

        AR_params.threshold = GL.get_threshold(edge_spread,signal_params.N);
        AR_params.gamma = GL.get_gamma(sparsities(s),AR_params.threshold,signal_params.N,vertice_std);
    
        f = zeros(trials,1);
        parfor i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,p);
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
        
            L = GL.AR_mean(y_noisy,AR_params);
    
            L_tmp = GL.threshold(L,AR_params.threshold);
            [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
        end
        F(s) = mean(f);
    end
end
