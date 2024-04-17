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
er_param = .1:.05:.95;
ba_param = 1:10;
rnd_param = .1:.1:1.2;
thresholds = 0:.0025:2;
gammas = logspace(-3,2,30);
%% N=20
signal_params = signals.create_default(20,'er');
signal_params.min_delay = 0;
signal_params.max_delay = 0;
[F_ER_20,S_ER_20,P_ER_20,R_ER_20,N_ER_20] = run_maximum_sims_dong(signal_params,er_param,trials,gammas,thresholds);
signal_params = signals.create_default(20,'gaussian');
signal_params.min_delay = 0;
signal_params.max_delay = 0;
[F_RND_20,S_RND_20,P_RND_20,R_RND_20,N_RND_20] = run_maximum_sims_dong(signal_params,rnd_param,trials,gammas,thresholds);
signal_params = signals.create_default(20,'pa');
signal_params.min_delay = 0;
signal_params.max_delay = 0;
[F_BA_20,S_BA_20,P_BA_20,R_BA_20,N_BA_20] = run_maximum_sims_dong(signal_params,ba_param,trials,gammas,thresholds);
%%
save("maximum_data_diff_density_dong_rect_FINAL_nodelay.mat")