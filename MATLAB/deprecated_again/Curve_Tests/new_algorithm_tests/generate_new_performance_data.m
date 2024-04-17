clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 64;
er_param = .1:.05:.95;
ba_param = 1:10;
rnd_param = .1:.1:1.2;
%% ER,N=40
% er_param = .1:.1:.9;
% signal_params = signals.create_default(40,'er');
% [F_ER_40,S_ER_40,T_ER_40] = run_performance_sims(signal_params,er_param,trials);
% er_param = .1:.05:.95;
%% RND,N=20
% signal_params = signals.create_default(20,'gaussian');
% [F_RND_20,S_RND_20,T_RND_20] = run_performance_sims(signal_params,rnd_param,trials);
%% ER,N=20
er_param = .1:.1:.9;
signal_params = signals.create_default(20,'er');
[F_ER_20,S_ER_20,T_ER_20] = run_performance_sims(signal_params,er_param,trials);
%% BA,N=20
% signal_params = signals.create_default(20,'pa');
% [F_BA_20,S_BA_20,T_BA_20] = run_performance_sims(signal_params,ba_param,trials);
%%
ba_param = 1:5;
%% RND,N=10
% signal_params = signals.create_default(10,'gaussian');
% [F_RND_10,S_RND_10,T_RND_10] = run_performance_sims(signal_params,rnd_param,trials);
% %% ER,N=10
% signal_params = signals.create_default(10,'er');
% [F_ER_10,S_ER_10,T_ER_10] = run_performance_sims(signal_params,er_param,trials);
%% BA,N=10
% signal_params = signals.create_default(10,'pa');
% [F_BA_10,S_BA_10,T_BA_10] = run_performance_sims(signal_params,ba_param,trials);
%%
save("large_new_performance_data_results_z_n20_test.mat")
%% 
function [F,sparsities,T] = run_performance_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);
    F = zeros(length(graph_params),1);
    sparsities = zeros(length(graph_params),1);

    T = zeros(length(graph_params),2);
    thresholds = 0:.0025:2;
    t_max = length(thresholds);

    for s = 1:length(graph_params)
        fprintf("%d %s",signal_params.N,signal_params.graph_type)
        s/length(graph_params)*100
        toc
        tic
    
        p = graph_params(s);
        sparsities(s) = graphs.get_sparsity(signal_params.graph_type,p,signal_params.N);
        sparsity = sparsities(s);

        a = -0.2724;
        b = 2.6540;
        g = a*exp(b.*sparsity);
        AR_params.gamma = 10^g;

        max_edges = signal_params.N*(signal_params.N-1)/2;
        num_edges = ceil(sparsity*max_edges);
        pred_edges = max(1,num_edges-1);

        f = zeros(trials,1);
        
        diff_t = zeros(trials,1);
        parfor i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,p);
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
            y_noisy = signals.z_score(y_noisy);

            L = GL.AR_mean(y_noisy,AR_params);

            % ft = zeros(t_max,1);
            % weight_count = zeros(t_max,1);
            % for tt = 1:t_max 
            %     L_tmp = GL.threshold(L,thresholds(tt));
            %     [~,~,ft(tt),~,~] = graphs.performance(L_0,L_tmp);
            %     weight_count(tt) = nnz(L_tmp-diag(diag(L_tmp)))/2;
            % end
            % indices = find(ft == max(ft));
            % diff_t(i) = thresholds(round(max(indices)));
            % diff_t(i) = weight_count(max(indices));

            weights = zeros(signal_params.N*(signal_params.N-1)/2,1);
            
            count = 1;
            for row = 1:(signal_params.N-1)
                for col = row+1:signal_params.N
                    weights(count) = L(row,col);
                    count = count+1;
                end
            end
            sorted_weights = sort(-weights,'descend');

            t = sorted_weights(pred_edges);
            
            if sparsity > .4 && sparsity < .6
                t = t*(1-.5*(sparsity-.4));
            elseif sparsity >= .6
                t = t*.9;
            end

            diff_t(i) = t;
            
            L_tmp = GL.threshold(L,t);
            [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);

        end
        F(s) = mean(f);
        T(s,1) = mean(diff_t);
        T(s,2) = pred_edges;
    end
end
