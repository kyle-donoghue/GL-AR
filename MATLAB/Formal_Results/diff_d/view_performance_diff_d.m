clc;clear;close all hidden;
addpath('../plot_functions');
%%
figure;
fig.fig = gcf;
%%
load('../../Formal_Data/diff_d/performance_data_diff_d_rect_FINAL.mat')
Q2= Q;
F_ER_202 = F_ER_20;
F_BA_202 = F_BA_20;
F_RND_202 = F_RND_20;
% load('../../Formal_Data/diff_d/performance_data_diff_d_rect_FINAL_40.mat');

fig.p1 = plot_perf(1,1,0:40,[F_ER_20(:,3);F_ER_20(1:end-1,3)],3);
fig.p2 = plot_perf(1,2,0:40,[F_BA_20(:,3);F_BA_20(1:end-1,3)],3);
fig.p3 = plot_perf(1,3,0:40,[F_RND_20(:,3);F_RND_20(1:end-1,3)],3);
%%
load('../../Formal_Data/diff_d/maximum_data_diff_d_dong_rect_FINAL.mat');
Q2= Q;
F_ER_surfs2 = F_ER_surfs;
F_BA_surfs2 = F_BA_surfs;
F_RND_surfs2 = F_RND_surfs;
load('../../Formal_Data/diff_d/maximum_data_diff_d_dong_rect_FINAL_40.mat');
Q = [Q2 Q];
F_ER_surfs2(:,:,end+1:end+20) = F_ER_surfs;
F_BA_surfs2(:,:,end+1:end+20) = F_BA_surfs;
F_RND_surfs2(:,:,end+1:end+20) = F_RND_surfs;
fig.p4 = plot_max(1,1,Q,F_ER_surfs2,3);
fig.p5 = plot_max(1,2,Q,F_BA_surfs2,3);
fig.p6 = plot_max(1,3,Q,F_RND_surfs2,3);
%%
subplot(1,3,1);
legend('GL-AR','GL-SigRep','GL-SigRep, E$[d]=10$','Interpreter','latex')
grid on;
STANDARDIZE_FIGURE(fig)
subplot(1,3,2);
legend('GL-AR','GL-SigRep','GL-SigRep, E$[d]=10$','Interpreter','latex')
grid on;
STANDARDIZE_FIGURE(fig)
subplot(1,3,3);
legend('GL-AR','GL-SigRep','GL-SigRep, E$[d]=10$','Interpreter','latex')
grid on;
STANDARDIZE_FIGURE(fig)