clc;clear;close all hidden;
load("gamma_threshold_sweep_with_sparsity_big_n10_ba.mat");
addpath(genpath('../../GL_classes/'));

%%
calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(signal_params.N*(signal_params.N-1)/2);
for i = 1:(length(sparsities)-2)
    figure(1);
    clf;hold on;
    [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds*signal_params.N);
    s = surf(X_axis,Y_axis,sparse_F(:,:,i)');
    s.EdgeColor = "none";
    view(2)
    ylim([0,10])
    [a,b] = GL.get_cliff_curve(sparsities(i));
    [c1,c2] = GL.plot_cliff_curve(gammas,a,b,signal_params.N);
    c1 = c1(c1<thresholds(end)*signal_params.N);
    c2 = c2(c2<thresholds(end)*signal_params.N);

    plot3(log10(gammas(1:length(c1))),c1,ones(length(c1),1),'LineWidth',2);
    plot3(log10(gammas(1:length(c2))),c2,ones(length(c2),1),'LineWidth',2);
    

    % edge_spread = graphs.get_edge_spread('er',sparsities(i),signal_params.N);
    edge_spread = graphs.get_edge_spread('pa',M(i),signal_params.N);

    t = GL.get_threshold(edge_spread,signal_params.N)*(signal_params.N);

    plot3(log10(gammas),t*ones(length(gammas),1),ones(length(gammas),1),'LineWidth',2,'Color','black');

    g = GL.get_gamma(sparsities(i),t/signal_params.N,signal_params.N);
    g = log10(g);

    plot3(g*ones(length(thresholds),1),thresholds*signal_params.N,ones(length(thresholds),1),'LineWidth',2,'Color','black');

    xlabel("Gamma")
    ylabel("Threshold")

    pause(1)
end
%%
