clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
%%
signal_params.N = 20;
signal_params.raw = 1;

signal_params.SNR = 2e0;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;

signal_params.intervals = signal_params.N;
signal_params.interval_length = 512;
signal_params.M = signal_params.intervals*signal_params.interval_length;

signal_params.order = 1;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .1;
signal_params.pole_mean = .5;

signal_params.min_delay = 10;
signal_params.max_delay = 50;


signal_params.minSeparation = 100;
signal_params.maxSeparation = 500;

trials = 1e2;
%% RND Generation
signal_params.graph_type = 'gaussian';
rnd_param = .05:.05:1.2;
variance_RND = zeros(length(rnd_param),trials);
for s = 1:length(rnd_param)
    s
    for i = 1:trials
        [L_0,A] = graphs.create(signal_params,rnd_param(s));
        variance_RND(s,i) = var(diag(L_0));
    end
end
variance_RND = mean(variance_RND,2,"omitnan");
std_RND = sqrt(variance_RND);
%% ER Generation
signal_params.graph_type = 'er';
er_param = .05:.05:.95;
variance_ER = zeros(length(er_param),trials);
for s = 1:length(er_param)
    s
    for i = 1:trials
        [L_0,A] = graphs.create(signal_params,er_param(s));
        variance_ER(s,i) = var(diag(L_0));
    end
end
variance_ER = mean(variance_ER,2,"omitnan");
std_ER = sqrt(variance_ER);
%% BA Generation
signal_params.graph_type = 'pa';
ba_param = 1:8;
variance_BA = zeros(length(ba_param),trials);
for s = 1:length(ba_param)
    s
    for i = 1:trials
        [L_0,A] = graphs.create(signal_params,ba_param(s));
        variance_BA(s,i) = var(diag(L_0));
    end
end
variance_BA = mean(variance_BA,2,"omitnan");
std_BA = sqrt(variance_BA);