clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('large_maximum_data_results_n20_3.mat');
load('large_maximum_data_results_n15_n10_3.mat');
%% N=10
figure(1);
subplot(1,3,1);
compare_threshold(F_ER_10,S_ER_10,10,er_param,'er','N10 ER')
subplot(1,3,2);
compare_threshold(F_BA_10,S_BA_10,10,1:8,'pa','N10 BA')
subplot(1,3,3);
compare_threshold(F_RND_10,S_RND_10,10,rnd_param,'gaussian','N10 RND')
%% N=15
figure(2);
subplot(1,3,1);
compare_threshold(F_ER_15,S_ER_15,15,er_param,'er','N15 ER')
subplot(1,3,2);
compare_threshold(F_BA_15,S_BA_15,15,1:11,'pa','N15 BA')
subplot(1,3,3);
compare_threshold(F_RND_15,S_RND_15,15,rnd_param,'gaussian','N15 RND')
%% N=20
figure(3);
subplot(1,3,1);
compare_threshold(F_ER_20,S_ER_20,20,er_param,'er','N20 ER')
subplot(1,3,2);
compare_threshold(F_BA_20,S_BA_20,20,ba_param,'pa','N20 BA')
subplot(1,3,3);
compare_threshold(F_RND_20,S_RND_20,20,rnd_param,'gaussian','N20 RND')
% %% N=40
% load('large_maximum_data_results_n40.mat');
% figure(3);
% subplot(1,3,1);
% compare_threshold(F_ER_40,S_ER_40,40,er_param,'er','N40 ER')
% subplot(1,3,2);
% compare_threshold(F_BA_40,S_BA_40,40,ba_param,'pa','N40 BA')
% subplot(1,3,3);
% compare_threshold(F_RND_40,S_RND_40,40,rnd_param,'gaussian','N40 RND')
%%
function compare_threshold(sparse_F,sparsities,N,graph_param,graph_type,tle)
    thresholds = 0:.0025:2;
    edge_spreads = graphs.get_edge_spread(graph_type,graph_param,N);
    xaxis = edge_spreads/N;
    max_thresholds_surf = zeros(length(sparsities),length(thresholds));
    max_thresholds = zeros(length(sparsities),1);
    z = zeros(length(sparsities),1);
    for i = 1:length(sparsities)
        maxes = max(sparse_F(:,:,i),[],1);
        max_thresholds_surf(i,:) = maxes;
        [z(i), ind] = max(sparse_F(:,:,i),[],"all");
        [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
        max_thresholds(i) = thresholds(col);
    end
    max_thresholds = max_thresholds*N;
    yaxis = thresholds*N;
    hold on;
    title(tle)
    xlabel('Normalized Edge Spread')
    ylabel('Normalized Threshold')
    zlabel('Fscore')
    [X_axis,Y_axis] = meshgrid(xaxis,yaxis);
    s = surf(X_axis,Y_axis,max_thresholds_surf');
    s.EdgeColor = "none";
    view(2)
    ylim([0 20])
    plot3(xaxis,max_thresholds,z)
    max_fit_thresholds = GL.get_threshold(edge_spreads,N)*N;
    plot3(xaxis,max_fit_thresholds,ones(size(xaxis)),'--')
end