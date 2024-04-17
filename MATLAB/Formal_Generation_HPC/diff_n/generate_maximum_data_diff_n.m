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
thresholds = 0:.0025:2;
gammas = logspace(-3,2,30);
N_list = [100:-4:52];
% N_list = 20;
%%
graph_param = rnd_param;
graph_type = 'gaussian';
%%
F_surfs = zeros(length(gammas),length(thresholds),length(N_list));
P_surfs = zeros(length(gammas),length(thresholds),length(N_list));
R_surfs = zeros(length(gammas),length(thresholds),length(N_list));
N_surfs = zeros(length(gammas),length(thresholds),length(N_list));
for i = 1:length(N_list)
    N = N_list(i);
    signal_params = signals.create_default(N,graph_type);
    signal_params.min_delay = 0;
    signal_params.max_delay = 4;
    [F_surfs(:,:,i),~,P_surfs(:,:,i),R_surfs(:,:,i),N_surfs(:,:,i)] = run_maximum_sims_dong(signal_params,graph_param,trials,gammas,thresholds);
    i/length(N_list)*100
    toc
    tic
end
%%
save("maximum_data_RND_dong_sine_FINAL_2delay_100.mat")
quit
