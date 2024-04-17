clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 32;
er_param = .2;
thresholds = 0:.0025:2;
gammas = logspace(-4,2,30);
snr_list = -10:1:20;
snr_list = 15;
snr_list = 10.^(snr_list/10);
N =20;
%%
F_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,'er');
    signal_params.minDelay = 0;
    signal_params.maxDelay = 0;
    signal_params.SNR=  snr_list(i);
    [F_surfs(:,:,i),~] = run_maximum_sims(signal_params,er_param,trials,gammas,thresholds);
    i/length(snr_list)*100
    toc
    tic
end
%%
% save("../../Formal_Data/diff_n/maximum_data_diff_n_PARAMETERS.mat")
