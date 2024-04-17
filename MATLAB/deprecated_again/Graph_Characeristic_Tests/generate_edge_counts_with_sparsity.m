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

trials = 1e3;
%% RND Generation
rnd_param = .05:.05:1.2;
edge_spread_RND = zeros(length(rnd_param),trials);
sparsity_RND = zeros(length(rnd_param),trials);
for s = 1:length(rnd_param)
    for i = 1:trials
        [L_0,A] = graphs.create(signal_params,'gaussian',0.75,rnd_param(s));
        num_edges = graphs.num_edges(L_0);
        sparsity_RND(s,i) = num_edges/(signal_params.N*(signal_params.N-1)/2);
        curr_max = max(diag(L_0));
        curr_avg = mean(diag(L_0));
        edge_spread_RND(s,i) = curr_avg/curr_max+curr_avg;
    end
end
edge_spread_RND = mean(edge_spread_RND,2,"omitnan");
sparsity_RND = mean(sparsity_RND,2);
%% ER Generation
er_param = .05:.05:.95;
edge_spread_ER = zeros(length(er_param),trials);
sparsity_ER = er_param';
for s = 1:length(er_param)
    for i = 1:trials
        [L_0,A] = graphs.create(signal_params,'er',er_param(s));
        curr_max = max(diag(L_0));
        curr_avg = mean(diag(L_0));
        edge_spread_ER(s,i) = curr_avg/curr_max+curr_avg;
    end
end
edge_spread_ER = mean(edge_spread_ER,2,"omitnan");
%% BA Generation
ba_param = 1:12;
edge_spread_BA = zeros(length(ba_param),trials);
sparsity_BA = zeros(length(ba_param),trials);
for s = 1:length(rnd_param)
    for i = 1:trials
        s
        i
        [L_0,A] = graphs.create(signal_params,'pa',ba_param(s));
        num_edges = graphs.num_edges(L_0);
        sparsity_BA(s,i) = num_edges/(signal_params.N*(signal_params.N-1)/2);
        curr_max = max(diag(L_0));
        curr_avg = mean(diag(L_0));
        edge_spread_BA(s,i) = curr_avg/curr_max+curr_avg;
    end
end
edge_spread_BA = mean(edge_spread_BA,2);
sparsity_BA = mean(sparsity_BA,2);