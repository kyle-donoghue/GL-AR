clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 128;
er_param = .1:.05:.95;
ba_param = 1:12;
rnd_param = .1:.1:1.2;
%% RND,N=20
signal_params = signals.create_default(20,'gaussian');
[F_RND_20,S_RND_20] = run_performance_sims(signal_params,rnd_param,trials);
%% ER,N=20
signal_params = signals.create_default(20,'er');
[F_ER_20,S_ER_20] = run_performance_sims(signal_params,er_param,trials);
%% BA,N=20
signal_params = signals.create_default(20,'pa');
[F_BA_20,S_BA_20] = run_performance_sims(signal_params,ba_param,trials);
%% ER,N=10
signal_params = signals.create_default(10,'er');
[F_ER_10,S_ER_10] = run_performance_sims(signal_params,er_param,trials);
%% BA,N=10
signal_params = signals.create_default(10,'pa');
ba_param = 1:8;
[F_BA_10,S_BA_10] = run_performance_sims(signal_params,ba_param,trials);
ba_param = 1:12;
%% RND,N=10
signal_params = signals.create_default(10,'gaussian');
[F_RND_10,S_RND_10] = run_performance_sims(signal_params,rnd_param,trials);
%% ER,N=15
signal_params = signals.create_default(15,'er');
[F_ER_15,S_ER_15] = run_performance_sims(signal_params,er_param,trials);
%% BA,N=15
signal_params = signals.create_default(15,'pa');
ba_param = 1:11;
[F_BA_15,S_BA_15] = run_performance_sims(signal_params,ba_param,trials);
ba_param = 1:12;
%% RND,N=15
signal_params = signals.create_default(15,'gaussian');
[F_RND_15,S_RND_15] = run_performance_sims(signal_params,rnd_param,trials);
%% ER,N=40
% signal_params = signals.create_default(40,'er');
% [F_ER_40,S_ER_40] = run_performance_sims(signal_params,er_param,trials);
% %% BA,N=40
% signal_params = signals.create_default(40,'pa');
% [F_BA_40,S_BA_40] = run_performance_sims(signal_params,ba_param,trials);
% %% RND,N=40
% signal_params = signals.create_default(40,'gaussian');
% [F_RND_40,S_RND_40] = run_performance_sims(signal_params,rnd_param,trials);
%%
save("large_performance_data_results__n20_n15_n10_6.mat")
%% 
function [F,sparsities] = run_performance_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);
    F = zeros(length(graph_params),1);
    sparsities = zeros(length(graph_params),1);
    for s = 1:length(graph_params)
        fprintf("%d %s",signal_params.N,signal_params.graph_type)
        s/length(graph_params)*100
        toc
        tic
    
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
