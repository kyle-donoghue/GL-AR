clc;clear;close all hidden;
load('large_maximum_data_results_n20_n10.mat');
load('large_maximum_data_results_n40.mat');
addpath(genpath('../../GL_classes/'));
%%
sparsities = S_RND_20;
sparse_F = F_RND_20;
N = 20;
graph_type = 'gaussian';
graph_param = rnd_param;
thresholds = 0:.0025:2;
gammas = logspace(-2,4,20);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
for i = 1:(length(sparsities))
    figure(1);
    clf;hold on;
    [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds*N);
    s = surf(X_axis,Y_axis,sparse_F(:,:,i)');
    s.EdgeColor = "none";
    view(2)
    ylim([0,20])
    [a,b] = GL.get_cliff_curve(sparsities(i));
    [c1,c2] = GL.plot_cliff_curve(gammas,a,b,N);
    c1 = c1(c1<thresholds(end)*N);
    c2 = c2(c2<thresholds(end)*N);

    plot3(log10(gammas(1:length(c1))),c1,ones(length(c1),1),'LineWidth',2);
    plot3(log10(gammas(1:length(c2))),c2,ones(length(c2),1),'LineWidth',2);
    

    edge_spread = graphs.get_edge_spread(graph_type,graph_param(i),N);

    t = GL.get_threshold(edge_spread,N)*(N);

    plot3(log10(gammas),t*ones(length(gammas),1),ones(length(gammas),1),'LineWidth',2,'Color','black');

    g = GL.get_gamma(sparsities(i),t/N,N);
    g = log10(g);

    plot3(g*ones(length(thresholds),1),thresholds*N,ones(length(thresholds),1),'LineWidth',2,'Color','black');

    xlabel("Gamma")
    ylabel("Threshold")

    pause(1)
end