clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 8;
er_param = .1:.1:.9;
ba_param = 1:8;
rnd_param = .1:.1:.8;
% %% ER,N=20
% signal_params = signals.create_default(20,'er');
% [W_ER_20,D_ER_20,T_ER_20] = run_weight_sims(signal_params,er_param,trials);
% %% RND,N=20
% signal_params = signals.create_default(20,'gaussian');
% [W_RND_20,D_RND_20,T_RND_20] = run_weight_sims(signal_params,rnd_param,trials);
%% ER,N=10
signal_params = signals.create_default(10,'er');
[W_ER_10,D_ER_10,T_ER_10] = run_weight_sims(signal_params,er_param,trials);
%% RND,N=10
signal_params = signals.create_default(10,'gaussian');
[W_RND_10,D_RND_10,T_RND_10] = run_weight_sims(signal_params,rnd_param,trials);
%% BA,N=10
signal_params = signals.create_default(10,'pa');
[W_BA_20,D_BA_20,T_BA_20] = run_weight_sims(signal_params,ba_param,trials);
%%
save("weights_n10.mat")
%% 
function [W,D,T] = run_weight_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);

    W = zeros(signal_params.N*(signal_params.N-1)/2,trials,length(graph_params));
    D = zeros(signal_params.N,trials,length(graph_params));
    T = zeros(trials,length(graph_params));
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
    
        max_t = zeros(trials,1);
        parfor i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,graph_params(s));
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
        
            L = GL.AR_mean(y_noisy,AR_params);
            f = zeros(t_max,1);
            weights = zeros(signal_params.N*(signal_params.N-1)/2,1);
            for t = 1:t_max 
                L_tmp = GL.threshold(L,thresholds(t));
                [~,~,f(t),~,~] = graphs.performance(L_0,L_tmp);
            end
            indices = find(f == max(f));
            max_t(i) = thresholds(round(mean(indices)));
            
            count = 1;
            for row = 1:(signal_params.N-1)
                for col = row+1:signal_params.N
                    weights(count) = L(row,col);
                    count = count+1;
                end
            end
            W(:,i,s) = weights;
            D(:,i,s) = diag(L);

        end
        T(:,s) = max_t;
    end
end
