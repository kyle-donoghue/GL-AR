clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 64;
er_param = .1:.05:.95;
ba_param = 1:10;
rnd_param = .1:.1:1.2;

%% N=20
signal_params = signals.create_default(20,'er');
AR_params = GL.create_default_params(signal_params);
[F_ER_20,S_ER_20,T_ER_20] = run_performance_sims(signal_params,AR_params,er_param,trials);
signal_params = signals.create_default(20,'gaussian');
AR_params = GL.create_default_params(signal_params);
[F_RND_20,S_RND_20,T_RND_20] = run_performance_sims(signal_params,AR_params,er_param,trials);
signal_params = signals.create_default(20,'pa');
AR_params = GL.create_default_params(signal_params);
[F_BA_20,S_BA_20,T_BA_20] = run_performance_sims(signal_params,AR_params,er_param,trials);

%%
save("../../Formal_Data/diff_density/performance_data_diff_density_PARAMETERS.mat")

