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
%%
trials = 500;
num_edges = zeros(trials,1);
% for i = 1:trials
%     [L_0,A] = graphs.create(signal_params,'gaussian',0.75,0.5);
%     [~,~,~,~,num_edges(i)] = graphs.performance(L_0,L_0);
% end
% graphs.plot(L_0,'L_0')
% mean(num_edges)/(signal_params.N*(signal_params.N-1)/2)
signal_params.graph_type = 'er';
graph_params = .05:.05:.95;
for s = 1:length(graph_params)
    for i = 1:trials
        [L_0,A] = graphs.create(signal_params,graph_params(s));
        [~,~,~,~,num_edges(i)] = graphs.performance(L_0,L_0);
        max_edges(i) = max(diag(L_0));
        avg_edges(i) = (max(diag(L_0))+min(diag(L_0)))/2;
    end
    avges(s) = mean(avg_edges);
    maxes(s) = mean(max_edges);
    sparsities(s) = mean(num_edges)/(signal_params.N*(signal_params.N-1)/2);
end
sparsities
maxes
avges