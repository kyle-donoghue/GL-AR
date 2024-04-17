clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 32;
er_param = .2;
thresholds = 0:.0025:2;
gammas = logspace(-4,2,30);
N_list = 100:-4:6;
N_list = 20;
%%
F_surfs = zeros(length(gammas),length(thresholds),length(N_list));
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,'er');
    signal_params.minDelay = 0;
    signal_params.maxDelay = 0;
    [F_surfs(:,:,i),~] = run_maximum_sims(signal_params,er_param,trials,gammas,thresholds);
    i/length(N_list)*100
    toc
    tic
end
%%
save("../../Formal_Data/sparsity_l1_trial/maximum_data_regular.mat")
