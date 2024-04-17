clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 8;
er_param = .1:.05:.95;
ba_param = 1:8;
rnd_param = .1:.1:1.2;
% %% ER,N=20
% signal_params = signals.create_default(20,'er');
% [T_ER_20] = run_threshold_sims(signal_params,er_param,trials);
% %% RND,N=20
% signal_params = signals.create_default(20,'gaussian');
% [T_RND_20] = run_threshold_sims(signal_params,rnd_param,trials);
% %% BA,N=20
% signal_params = signals.create_default(20,'pa');
% [T_BA_20] = run_threshold_sims(signal_params,ba_param,trials);
%% ER,N=10
signal_params = signals.create_default(10,'er');
[T_ER_10] = run_threshold_sims(signal_params,er_param,trials);
%% RND,N=10
signal_params = signals.create_default(10,'gaussian');
[T_RND_10] = run_threshold_sims(signal_params,rnd_param,trials);
%% BA,N=10
signal_params = signals.create_default(10,'pa');
[T_BA_10] = run_threshold_sims(signal_params,ba_param,trials);
%%
save("thresholds_n10.mat")
%% 
function [T] = run_threshold_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);

    T = zeros(length(graph_params),1);

    thresholds = 0:.0025:2;
    t_max = length(thresholds);

    for s = 1:length(graph_params)
        fprintf("%d %s",signal_params.N,signal_params.graph_type)
        s/length(graph_params)*100
        toc
        tic

        sparsity = graphs.get_sparsity(signal_params.graph_type,graph_params(s));

        a = -67.8256;
        b = 0.0413;
        c = 69.3318;
        g = c+a*exp(b.*sparsity.^2);

        AR_params.gamma = 10^g;
    
        num_edges = ceil(sparsity*signal_params.N*(signal_params.N-1)/2);

        max_t = zeros(trials,1);
        for i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,graph_params(s));
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
        
            L = GL.AR_mean(y_noisy,AR_params);
            f = zeros(t_max,1);

            weights = zeros(signal_params.N*(signal_params.N-1)/2,1);
            
            count = 1;
            for row = 1:(signal_params.N-1)
                for col = row+1:signal_params.N
                    weights(count) = L(row,col);
                    count = count+1;
                end
            end
            sorted_weights = sort(-weights,'descend')*signal_params.N;
            
            for nn = (1+num_edges):length(sorted_weights)
                diff = abs(sorted_weights(nn-1)-sorted_weights(nn))/sorted_weights(nn-1);
                if sorted_weights(nn-1) < 1e-3
                    index = nn-1;
                    break
                end
                if diff > .1
                    index = nn-1;
                    break;
                end
            end

            max_t(i) = sorted_weights(index)

        end
        T(s) = mean(max_t);
    end
end
