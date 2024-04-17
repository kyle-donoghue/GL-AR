clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load("large_maximum_data_results_n20_3.mat");
load("large_maximum_data_results_n15_n10_3.mat");
%%
figure;hold on;
%%
signal_params.N = 20;
sparsities = S_ER_20;
sparse_F = F_ER_20;
graph_type = 'er';
graph_param = er_param;
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
plot(xaxis,max_thresholds,'LineWidth',2)
xlabel('Number of Edges')
ylabel('Threshold*N that Produces Fmax')
%%
sparsities = S_BA_20;
sparse_F = F_BA_20;
graph_type = 'pa';
graph_param = ba_param;
gammas = logspace(-2,4,60);
thresholds = 0:.0025:2;
spread_edges = graphs.get_edge_spread2(graph_type,ba_param,signal_params.N)';
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
plot(xaxis,max_thresholds,'LineWidth',2)
xlabel('Edge Spread')
ylabel('Threshold*N that Produces Fmax')
%%
sparsities = S_RND_20;
sparse_F = F_RND_20;
graph_type = 'gaussian';
graph_param = rnd_param;
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
plot(xaxis,max_thresholds,'LineWidth',2)
%%
edge_spread2 = .01:.01:300;
plot(edge_spread2,GL.get_threshold(edge_spread2,20)*20,'--','LineWidth',1.5)
%%
legend('ER','BA','RND','Fit')
xlabel('Edge Spread, $$\lambda$$','Interpreter','latex','FontSize',14)
ylabel('Optimal $$\tilde{\delta}$$','Interpreter','latex','FontSize',14)
% title('Normalized Threshold that Produces Maximum F-Score','Interpreter','latex','FontSize',14)
title('N=20','Interpreter','latex','FontSize',14)
%%
figure; hold on;
sparsities = S_RND_20;
sparse_F = F_RND_20;
graph_type = 'gaussian';
graph_param = rnd_param;
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
plot(xaxis,max_thresholds)
%%
signal_params.N = 10;
sparsities = S_RND_10;
sparse_F = F_RND_10;
graph_type = 'gaussian';
graph_param = rnd_param;
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
plot(xaxis,max_thresholds)
%%
signal_params.N = 15;
sparsities = S_RND_15;
sparse_F = F_RND_15;
graph_type = 'gaussian';
graph_param = rnd_param;
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
plot(xaxis,max_thresholds)
%%
edge_spread2 = .01:.01:300;
plot(edge_spread2,GL.get_threshold(edge_spread2,20)*20,'--')
edge_spread2 = .01:.01:300;
plot(edge_spread2,GL.get_threshold(edge_spread2,15)*15,'--')
edge_spread2 = .01:.01:300;
plot(edge_spread2,GL.get_threshold(edge_spread2,10)*10,'--')
%%
legend('rnd20','rnd10','rnd15','fit20','fit15','fit10')
xlabel('Edge Spread')
ylabel('Threshold*N that Produces Fmax')
%%
figure;hold on;
%%
signal_params.N = 15;
sparsities = S_ER_15;
sparse_F = F_ER_15;
graph_type = 'er';
graph_param = er_param;
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
plot(xaxis,max_thresholds,'LineWidth',2)
xlabel('Number of Edges')
ylabel('Threshold*N that Produces Fmax')
%%
sparsities = S_BA_15;
sparse_F = F_BA_15;
graph_type = 'pa';
graph_param = 1:11;
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
plot(xaxis,max_thresholds,'LineWidth',2)
xlabel('Edge Spread')
ylabel('Threshold*N that Produces Fmax')
%%
sparsities = S_RND_15;
sparse_F = F_RND_15;
graph_type = 'gaussian';
graph_param = rnd_param;
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
plot(xaxis,max_thresholds,'LineWidth',2)
%%
edge_spread2 = .01:.01:300;
plot(edge_spread2,GL.get_threshold(edge_spread2,15)*15,'--','LineWidth',1.5)
%%
legend('ER','BA','RND','Fit')
xlabel('Edge Spread, $$\lambda$$','Interpreter','latex','FontSize',14)
ylabel('Optimal $$\tilde{\delta}$$','Interpreter','latex','FontSize',14)
title('N=15','Interpreter','latex','FontSize',14)
%%
figure;hold on;
%%
signal_params.N = 10;
sparsities = S_ER_10;
sparse_F = F_ER_10;
graph_type = 'er';
graph_param = er_param;
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
plot(xaxis,max_thresholds,'LineWidth',2)
xlabel('Number of Edges')
ylabel('Threshold*N that Produces Fmax')
title('10')
%%
sparsities = S_BA_10;
sparse_F = F_BA_10;
graph_type = 'pa';
graph_param = 1:8;
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
plot(xaxis,max_thresholds,'LineWidth',2)
xlabel('Edge Spread')
ylabel('Threshold*N that Produces Fmax')
%%
sparsities = S_RND_10;
sparse_F = F_RND_10;
graph_type = 'gaussian';
graph_param = rnd_param;
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
plot(xaxis,max_thresholds,'LineWidth',2)
%%
edge_spread2 = .01:.01:300;
plot(edge_spread2,GL.get_threshold(edge_spread2,10)*10,'--','LineWidth',1.5)
%%
legend('ER','BA','RND','Fit')
xlabel('Edge Spread, $$\lambda$$','Interpreter','latex','FontSize',14)
ylabel('Optimal $$\tilde{\delta}$$','Interpreter','latex','FontSize',14)
title('N=10','Interpreter','latex','FontSize',14)