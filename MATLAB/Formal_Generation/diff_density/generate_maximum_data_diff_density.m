clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
addpath('../sim_functions/')
%%
trials = 32;
er_param = .1:.05:.95;
er_param = .2;
ba_param = 1:10;
rnd_param = .1:.1:1.2;

thresholds = 0:.0025:2;
% gammas = logspace(-3,2,30);
gammas = logspace(-17.5,-17,15);
%% N=20
signal_params = signals.create_default(20,'er');
% [F_ER_20,S_ER_20] = run_maximum_sims(signal_params,er_param,trials,gammas, thresholds);
[F_ER_20,S_ER_20] = run_maximum_L1_sims(signal_params,er_param,trials,gammas, thresholds);
% signal_params = signals.create_default(20,'gaussian');
% AR_params = GL.create_default_params(signal_params);
% [F_RND_20,S_RND_20] = run_maximum_sims(signal_params,AR_params,er_param,trials);
% signal_params = signals.create_default(20,'pa');
% AR_params = GL.create_default_params(signal_params);
% [F_BA_20,S_BA_20] = run_maximum_sims(signal_params,AR_params,er_param,trials);
%%
save("../../Formal_Data/diff_density/maximum_data_diff_density_L1_default2.mat")