clc;clear;close all hidden;
load("gamma_threshold_sweep_with_sparsity_big.mat");
%%
max_gammas = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_gammas(i) = gammas(row);
end
%%
figure;hold on;
plot(sparsities,log10(max_gammas))
mean(log10(max_gammas))
%%
% fit_data = fit(sparsities'-sparsities(1),max_gammas-1/(2*signal_params.N),'exp1');
% y_hat = fit_data.a*exp(fit_data.b*(sparsities-sparsities(1)))+1/(2*signal_params.N);
% plot(sparsities,y_hat)