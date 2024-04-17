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
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds(14:end) = .75*max_thresholds(13);
% max_thresholds = max_thresholds*signal_params.N;
%%
% num_edges = sparsities*signal_params.N;
num_edges = sparsities*(signal_params.N*(signal_params.N-1)/2);
%%
figure;hold on;
plot(num_edges,max_thresholds)
%%
weights = ones(length(num_edges),1);
% weights([8 13]) = 10;
weights([8]) = 5;
fit_data = fit(num_edges'-num_edges(1),max_thresholds-max_thresholds(end),'exp1','Weights',weights,'StartPoint',[0,0])
y_hat = fit_data.a*exp(fit_data.b*(num_edges-num_edges(1)))+max_thresholds(end);
plot(num_edges,y_hat)
xlabel('Number of Edges')
ylabel('Threshold that Produces Fmax')