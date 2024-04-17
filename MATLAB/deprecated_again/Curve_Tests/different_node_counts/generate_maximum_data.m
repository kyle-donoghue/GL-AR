clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
trials = 32;
er_param = .1:.05:.95;
ba_param = 1:12;
rnd_param = .1:.1:1.2;
% %% ER,N=20
% signal_params = signals.create_default(20,'er');
% [F_ER_20,S_ER_20] = run_maximum_sims(signal_params,er_param,trials);
% %% BA,N=20
% signal_params = signals.create_default(20,'pa');
% [F_BA_20,S_BA_20] = run_maximum_sims(signal_params,ba_param,trials);
% %% RND,N=20
% signal_params = signals.create_default(20,'gaussian');
% [F_RND_20,S_RND_20] = run_maximum_sims(signal_params,rnd_param,trials);
% %% ER,N=15
% signal_params = signals.create_default(15,'er');
% [F_ER_15,S_ER_15] = run_maximum_sims(signal_params,er_param,trials);
% %% BA,N=15
% signal_params = signals.create_default(15,'pa');
% ba_param = 1:11;
% [F_BA_15,S_BA_15] = run_maximum_sims(signal_params,ba_param,trials);
% ba_param = 1:12;
% %% RND,N=15
% signal_params = signals.create_default(15,'gaussian');
% [F_RND_15,S_RND_15] = run_maximum_sims(signal_params,rnd_param,trials);
% %% ER,N=10
% signal_params = signals.create_default(10,'er');
% [F_ER_10,S_ER_10] = run_maximum_sims(signal_params,er_param,trials);
% %% BA,N=10
% signal_params = signals.create_default(10,'pa');
% ba_param = 1:8;
% [F_BA_10,S_BA_10] = run_maximum_sims(signal_params,ba_param,trials);
% ba_param = 1:12;
% %% RND,N=10
% signal_params = signals.create_default(10,'gaussian');
% [F_RND_10,S_RND_10] = run_maximum_sims(signal_params,rnd_param,trials);
% %% ER,N=40
% signal_params = signals.create_default(40,'er');
% [F_ER_40,S_ER_40] = run_maximum_sims(signal_params,er_param,trials);
% %% BA,N=40
% signal_params = signals.create_default(40,'pa');
% [F_BA_40,S_BA_40] = run_maximum_sims(signal_params,ba_param,trials);
% %% RND,N=40
% signal_params = signals.create_default(40,'gaussian');
% [F_RND_40,S_RND_40] = run_maximum_sims(signal_params,rnd_param,trials);
%% ER,N=100
signal_params = signals.create_default(100,'er');
[F_ER_100,S_ER_100] = run_maximum_sims(signal_params,.2,trials);
%%
save("large_maximum_data_results_n100.mat")
%% 
function [sparse_F,sparsities] = run_maximum_sims(signal_params,graph_params,trials)
    AR_params = GL.create_default_params(signal_params);
    
    thresholds = 0:.0025:20;
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
end
