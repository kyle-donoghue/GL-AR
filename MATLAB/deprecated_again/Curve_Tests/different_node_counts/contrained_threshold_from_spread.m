clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load("large_maximum_data_results_n15_n10_3.mat");
load("large_maximum_data_results_n20_3.mat");
%%

%%
% sparsities = S_ER_20;
% sparse_F = F_ER_20;
% graph_type = 'er';
% graph_param = er_param;

sparsities = S_RND_20;
sparse_F = F_RND_20;
graph_type = 'gaussian';
graph_param = rnd_param;

% sparsities = S_BA_20;
% sparse_F = F_BA_20;
% graph_type = 'pa';
% graph_param = ba_param;

% signal_params.N = 15;
% sparsities = S_RND_15;
% sparse_F = F_RND_15;
% graph_type = 'gaussian';
% graph_param = rnd_param;

% signal_params.N = 10;
% sparsities = S_RND_10;
% sparse_F = F_RND_10;
% graph_type = 'gaussian';
% graph_param = rnd_param;

gammas = logspace(-2,4,60);
thresholds = 0:.0025:2;

spread_edges = graphs.get_edge_spread2(graph_type,graph_param,signal_params.N)';
xaxis = spread_edges;
max_thresholds = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    v_std = graphs.get_vertice_std(graph_type,graph_param(i),signal_params.N);
    [a,b] = GL.get_cliff_curve(sparsities(i));
    [c1,~] = GL.plot_cliff_curve(gammas,a,b,signal_params.N,v_std);
    % c1 = c1(c1<thresholds(end)*signal_params.N);

    [X,Y] = meshgrid(log10(gammas),thresholds*signal_params.N);
    cliff_curve_z1 = interp2(X,Y,sparse_F(:,:,i)',log10(gammas),c1);
    
    [~, ind] = max(cliff_curve_z1,[],"all");
    max_thresholds(i) = c1(ind);
end
%%
figure;hold on;
plot(xaxis,max_thresholds)
ridge_offset = .9;
max_thresholds2 = max_thresholds - ridge_offset;
%%
weights = ones(length(xaxis),1);
% weights([8 13]) = 10;
% weights([8]) = 5;
o = fitoptions('exp2');
o.Lower = [0 -Inf 0 -Inf];
o.Upper = [Inf 0 Inf 0];
fit_data = fit(xaxis,max_thresholds-ridge_offset,'exp2',o)
y_hat = fit_data.a*exp(fit_data.b*(xaxis))+fit_data.c*exp(fit_data.d*(xaxis))+ridge_offset;
plot(xaxis,y_hat,'--')
xlabel('Number of Edges')
ylabel('Threshold*N that Produces Fmax')