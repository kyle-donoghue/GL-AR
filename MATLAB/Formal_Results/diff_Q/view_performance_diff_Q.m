clc;clear;close all hidden;
addpath('../plot_functions');
%%
figure;
fig.fig = gcf;
%%
load('../../Formal_Data/diff_Q/performance_data_diff_Q_rect_FINAL.mat')
fig.p1 = plot_perf(1,1,Q,F_ER_20(:,3),3);
fig.p2 = plot_perf(1,2,Q,F_BA_20(:,3),3);
fig.p3 = plot_perf(1,3,Q,F_RND_20(:,3),3);
%%
load('../../Formal_Data/diff_Q/maximum_data_diff_Q_dong_rect_FINAL_nodelay.mat');
fig.p4 = plot_max(1,1,Q,F_ER_surfs,3);
fig.p5 = plot_max(1,2,Q,F_BA_surfs,3);
fig.p6 = plot_max(1,3,Q,F_RND_surfs,3);
%%
load('../../Formal_Data/diff_Q/maximum_data_diff_Q_dong_rect_FINAL_10delay.mat');
fig.p4 = plot_max(1,1,Q,F_ER_surfs,3);
fig.p5 = plot_max(1,2,Q,F_BA_surfs,3);
fig.p6 = plot_max(1,3,Q,F_RND_surfs,3);
%%
subplot(1,3,1);
legend('GL-AR, E$[d]=10$','GL-SigRep, E$[d]=0$','GL-SigRep, E$[d]=10$','Interpreter','latex')
grid on;
STANDARDIZE_FIGURE(fig)
subplot(1,3,2);
legend('GL-AR, E$[d]=10$','GL-SigRep, E$[d]=0$','GL-SigRep, E$[d]=10$','Interpreter','latex')
grid on;
STANDARDIZE_FIGURE(fig)
subplot(1,3,3);
legend('GL-AR, E$[d]=10$','GL-SigRep, E$[d]=0$','GL-SigRep, E$[d]=10$','Interpreter','latex')
grid on;
STANDARDIZE_FIGURE(fig)