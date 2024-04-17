clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 64;
er_param = .2;
%% ER=.2
N_list = 50:-4:6;

% N_list = 20;

F_list = zeros(length(N_list),1);
T_list = zeros(length(N_list),2);
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,'er');
    signal_params.SNR = 200;
    % signal_params.directed = 0;
    [F_list(i),~,T_list(i,:)] = run_performance_sims(signal_params,er_param,trials);
    i/length(N_list)*100
    toc
    tic
end
%%
save("large_new_performance_data_results_diff_n_er_z_SNR200_AR100_quadprog.mat")
%% 
function [F,sparsities,T] = run_performance_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);
    AR_params.P = 100;

    F = zeros(length(graph_params),1);
    sparsities = zeros(length(graph_params),1);
    for s = 1:length(graph_params)
        fprintf("%d %s",signal_params.N,signal_params.graph_type)
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
        t_min = .1/signal_params.N;

        f = zeros(trials,1);
        diff_t = zeros(trials,1);
        parfor i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,p);
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
            y_noisy = signals.z_score(y_noisy);

            L = GL.AR_mean(y_noisy,AR_params);
    
            weights = zeros(signal_params.N*(signal_params.N-1)/2,1);
            
            count = 1;
            for row = 1:(signal_params.N-1)
                for col = row+1:signal_params.N
                    weights(count) = L(row,col);
                    count = count+1;
                end
            end
            sorted_weights = sort(-weights,'descend');
            % neg_index = find(sorted_weights < 0,1, 'first');
            % if isempty(neg_index)
            %     t = sorted_weights(pred_edges);
            % else
            %     t = t_min;
            % end
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
