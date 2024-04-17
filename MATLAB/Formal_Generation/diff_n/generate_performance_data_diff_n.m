clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 32;
er_param = .2;
rnd_param = .4;
% ba_param = .4;
% N_list = 100:-4:6;
N_list = 100:-2:10;
N_list = 20;
%%
graph_type = 'er';
graph_param = er_param;
%% Simulation
F_list = zeros(length(N_list),5);
T_list = zeros(length(N_list),2);
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,graph_type);
    AR_params = GL.create_default_params(signal_params);
    % AR_params.P = 75;
    % signal_params.SNR = 1000;
    signal_params.minDelay = 0;
    signal_params.maxDelay = 0;
    [F_list(i,:),~,T_list(i,:)] = run_performance_sims(signal_params,AR_params,graph_param,trials);
    i/length(N_list)*100
    toc
    tic
end
%%
save("../../Formal_Data/diff_n/performance_data_diff_n_RND_sine_FIyuNAL.mat")