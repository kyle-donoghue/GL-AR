clc;clear;close all hidden;
load("gamma_threshold_sweep_with_sparsity_big.mat");
%%
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
max_thresholds(1) = .5;
% max_thresholds(14:end) = 1/(2*signal_params.N);
max_thresholds(14:end) = .75/signal_params.N;
%%
figure;hold on;
plot(sparsities,max_thresholds)
%%
weights = ones(length(sparsities),1);
% weights([8 13]) = 10;
weights([8]) = 5;
fit_data = fit(sparsities'-sparsities(1),max_thresholds-max_thresholds(14),'exp1','Weights',weights)
y_hat = fit_data.a*exp(fit_data.b*(sparsities-sparsities(1)))+max_thresholds(14);
plot(sparsities,y_hat)