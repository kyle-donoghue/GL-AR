clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 64;
er_param = .2;
N_list = 50:-4:6;

%% Simulation
F_list = zeros(length(N_list),1);
T_list = zeros(length(N_list),2);
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,'er');
    AR_params = GL.create_default_params(signal_params);
    [F_list(i),~,T_list(i,:)] = run_performance_sims(signal_params,AR_params,er_param,trials);
    i/length(N_list)*100
    toc
    tic
end
%%
save("../../Formal_Data/diff_n/performance_data_diff_n_PARAMETERS.mat")