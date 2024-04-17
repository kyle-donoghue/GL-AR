clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
cores=32
poolobj = gcp('nocreate'); % If pool, delete
if ~isempty(poolobj)
    if poolobj.NumWorkers ~= cores
        delete(poolobj);
        poolobj = parpool('local',cores);
    end
else
    poolobj = parpool('local',cores);
end
fprintf('Number of workers: %g\n', poolobj.NumWorkers);
tic
%%
trials = 32;
er_param = .2;
rnd_param = .4;
ba_param = 2;
thresholds = 0:.0025:2;
gammas = logspace(-3,2,30);
snr_list = -10:1:20;
snr_list = 10.^(snr_list/10);
N = 20;
%%
graph_param = er_param;
graph_type = 'er';
F_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
P_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
R_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
N_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,graph_type);
    signal_params.min_delay = 5;
    signal_params.max_delay = 15;
    signal_params.SNR = snr_list(i);
    [F_surfs(:,:,i),~,P_surfs(:,:,i),R_surfs(:,:,i),N_surfs(:,:,i)] = run_maximum_sims_dong(signal_params,graph_param,trials,gammas,thresholds);
    i/length(snr_list)*100
    toc
    tic
end
F_ER_surfs = F_surfs;
P_ER_surfs = P_surfs;
R_ER_surfs = R_surfs;
N_ER_surfs = N_surfs;
%%
graph_param = ba_param;
graph_type = 'pa';
F_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
P_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
R_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
N_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,graph_type);
    signal_params.min_delay = 5;
    signal_params.max_delay = 15;
    signal_params.SNR = snr_list(i);
    [F_surfs(:,:,i),~,P_surfs(:,:,i),R_surfs(:,:,i),N_surfs(:,:,i)] = run_maximum_sims_dong(signal_params,graph_param,trials,gammas,thresholds);
    i/length(snr_list)*100
    toc
    tic
end
F_BA_surfs = F_surfs;
P_BA_surfs = P_surfs;
R_BA_surfs = R_surfs;
N_BA_surfs = N_surfs;
%%
graph_param = rnd_param;
graph_type = 'gaussian';
F_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
P_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
R_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
N_surfs = zeros(length(gammas),length(thresholds),length(snr_list));
for i = 1:length(snr_list)
    signal_params = signals.create_default(N,graph_type);
    signal_params.min_delay = 5;
    signal_params.max_delay = 15;
    signal_params.SNR = snr_list(i);
    [F_surfs(:,:,i),~,P_surfs(:,:,i),R_surfs(:,:,i),N_surfs(:,:,i)] = run_maximum_sims_dong(signal_params,graph_param,trials,gammas,thresholds);
    i/length(snr_list)*100
    toc
    tic
end
F_RND_surfs = F_surfs;
P_RND_surfs = P_surfs;
R_RND_surfs = R_surfs;
N_RND_surfs = N_surfs;
%%
save("maximum_data_diff_snr_dong_sine_FINAL_10delay_rect.mat")
quit
