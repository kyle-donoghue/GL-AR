clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 64;
er_param = .2;
rnd_param = .4;
ba_param = 2;
snr_list = -10:1:20;
snr_list = 10.^(snr_list/10);
N =20;
%%
graph_type = 'er';
graph_param = er_param;
F_list = zeros(length(snr_list),5);
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,graph_type);
    AR_params = GL.create_default_params(signal_params);
    signal_params.SNR = snr_list(i);
    % AR_params.P = 100;
    [F_list(i,:)] = run_performance_sims(signal_params,AR_params,graph_param,trials);
    i/length(snr_list)*100
    toc
    tic
end
F_ER_20 = F_list;
%%
graph_type = 'pa';
graph_param = ba_param;
F_list = zeros(length(snr_list),5);
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,graph_type);
    AR_params = GL.create_default_params(signal_params);
    signal_params.SNR = snr_list(i);
    [F_list(i,:)] = run_performance_sims(signal_params,AR_params,graph_param,trials);
    i/length(snr_list)*100
    toc
    tic
end
F_BA_20 = F_list;
%%
graph_type = 'gaussian';
graph_param = rnd_param;
F_list = zeros(length(snr_list),5);
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,graph_type);
    AR_params = GL.create_default_params(signal_params);
    signal_params.SNR = snr_list(i);
    [F_list(i,:)] = run_performance_sims(signal_params,AR_params,graph_param,trials);
    i/length(snr_list)*100
    toc
    tic
end
F_RND_20 = F_list;
%%
save("../../Formal_Data/diff_snr/performance_data_diff_snr_sine_FINAL_rect.mat")