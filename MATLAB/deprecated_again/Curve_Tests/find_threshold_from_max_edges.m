clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
load("gamma_threshold_sweep_with_sparsity_big.mat");
%%
% min_edges = zeros(length(sparsities),500);
% max_edges = zeros(length(sparsities),500);
% avg_edges = zeros(length(sparsities),500);
% for i = 1:length(sparsities)
%     for t = 1:500
%         [L_0,A] = graphs.create(signal_params,'er',sparsities(i));
%         maxes = sort(diag(L_0));
%         % max_edges(i,t) = max(diag(L_0));
%         min_edges(i,t) = min(diag(L_0));
%         max_edges(i,t) = sum(maxes((end-floor(length(maxes)/5)+1):end));
%         max_edges(i,t) = maxes(end);
%         % avg_edges(i,t) = (max(diag(L_0))+min(diag(L_0)))/2;
%         avg_edges(i,t) = mean(diag(L_0));
%     end
% end
% min_edges = mean(min_edges,2);
% max_edges = mean(max_edges,2);
% avg_edges = mean(avg_edges,2);
% sqrt_edges = log(avg_edges);
% sum_edges = (max_edges+floor(length(maxes)/2)*avg_edges)/(floor(length(maxes)/2)+1);
% spread_edges = (avg_edges)./max_edges+avg_edges;
spread_edges = graphs.get_edge_spread('er',sparsities,signal_params.N)';
xaxis = spread_edges/signal_params.N;
%%
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_thresholds(i) = thresholds(col);
end
% max_thresholds(1) = .5;
% % max_thresholds(14:end) = 1/(2*signal_params.N);
% % max_thresholds(14:end) = .75/signal_params.N;
% max_thresholds(14:end) = .75/signal_params.N;
max_thresholds = max_thresholds*signal_params.N;
%%
figure;hold on;
plot(xaxis,max_thresholds)
%%
weights = ones(length(xaxis),1);
% weights([8 13]) = 10;
% weights([8]) = 5;
fit_data = fit(xaxis,max_thresholds-max_thresholds(end),'exp1','Weights',weights,'StartPoint',[0,0])
y_hat = fit_data.a*exp(fit_data.b*(xaxis))+max_thresholds(end);
plot(xaxis,y_hat)
xlabel('Number of Edges')
ylabel('Threshold*N that Produces Fmax')