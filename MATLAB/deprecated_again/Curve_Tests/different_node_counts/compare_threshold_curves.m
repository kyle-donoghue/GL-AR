clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Original ER Plot
load('gamma_threshold_sweep_with_sparsity_big.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
max_thresholds(1) = .5;
max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure;hold on;
title('N20 ER')
plot(xaxis,max_thresholds)
% fit_data = fit(xaxis',max_thresholds-max_thresholds(end),'exp1','StartPoint',[0,0])
% y_hat = fit_data.a*exp(fit_data.b*(xaxis))+.75;
% plot(xaxis,y_hat,'-*')
%% Current Fit
xaxis2 = xaxis;
max_fit_thresholds = GL.get_threshold(edge_spreads,signal_params.N)*signal_params.N;
plot(xaxis2,max_fit_thresholds,'--')
%% Small ER Plot
load('gamma_threshold_sweep_with_sparsity_big_n10.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure; hold on;
title('N10 ER')
plot(xaxis2,max_fit_thresholds,'--')
plot(xaxis,max_thresholds)
%% Big BA Plot
load('gamma_threshold_sweep_with_sparsity_small_BA.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure; hold on;
title('N20 BA')
plot(xaxis2,max_fit_thresholds,'--')
plot(xaxis,max_thresholds)
%% Big RND Plot
load('gamma_threshold_sweep_with_sparsity_small_RND2.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure; hold on;
title('N20 RND')
plot(xaxis2,max_fit_thresholds,'--')
plot(xaxis,max_thresholds)
%% Small BA Plot
load('gamma_threshold_sweep_with_sparsity_big_n10_ba.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure; hold on;
title('N10 BA')
plot(xaxis2,max_fit_thresholds,'--')
plot(xaxis,max_thresholds)
%% Small RND Plot
load('gamma_threshold_sweep_with_sparsity_big_n10_rnd.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure; hold on;
title('N10 RND')
plot(xaxis2,max_fit_thresholds,'--')
plot(xaxis,max_thresholds)
%% Very Big ER Plot
load('gamma_threshold_sweep_with_sparsity_big_n40.mat');
edge_spreads = graphs.get_edge_spread('er',sparsities,signal_params.N);
xaxis = edge_spreads/signal_params.N;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
figure; hold on;
title('N40 ER')
plot(xaxis2,max_fit_thresholds,'--')
plot(xaxis,max_thresholds)