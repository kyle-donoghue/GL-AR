clc;clear;close all hidden;
addpath('../plot_functions');
%%
figure;
fig.fig = gcf;
load('../../Formal_Data/diff_density/performance_data_diff_density_FINAL.mat')
%% N = 20
plot_perf(1,1,S_ER_20,F_ER_20(:,3),3);
plot_perf(1,2,S_BA_20,F_BA_20(:,3),3);
plot_perf(1,3,S_RND_20(2:end),F_RND_20(2:end,3),3);
%%
%%
load('../../Formal_Data/diff_density/maximum_data_diff_density_dong_rect_FINAL_nodelay.mat');
fig.p4 = plot_max(1,1,S_ER_20,F_ER_20,3);
fig.p5 = plot_max(1,2,S_BA_20,F_BA_20,3);
fig.p6 = plot_max(1,3,S_RND_20(2:end),F_RND_20(:,:,2:end),3);
%%
load('../../Formal_Data/diff_density/initial_results.mat');
F_ER_202 = F_ER_20;
load('../../Formal_Data/diff_density/maximum_data_diff_density_dong_rect_FINAL_10delay.mat');
F_ER_20 = F_ER_202;
fig.p4 = plot_max(1,1,S_ER_20,F_ER_20,3);
fig.p5 = plot_max(1,2,S_BA_20,F_BA_20,3);
fig.p6 = plot_max(1,3,S_RND_20(2:end),F_RND_20(:,:,2:end),3);
%% N = 20
% plot_max(2,1,S_ER_20,F_ER_20);
% plot_max(1,2,(1:10)*2/20,F_BA_20(:,:,1:10));
% plot_max(1,3,S_RND_20,F_RND_20);
subplot(1,3,1);
% plot(S_RND_20,graphs.minFscore(S_RND_20),'--')
legend('GL-AR, E$[d]=10$','GL-SigRep, E$[d]=0$','GL-SigRep, E$[d]=10$','Interpreter','latex','Location','southeast')
grid on;
STANDARDIZE_FIGURE(fig)
subplot(1,3,2);
legend('GL-AR, E$[d]=10$','GL-SigRep, E$[d]=0$','GL-SigRep, E$[d]=10$','Interpreter','latex','Location','southeast')
grid on;
STANDARDIZE_FIGURE(fig)
subplot(1,3,3);
legend('GL-AR, E$[d]=10$','GL-SigRep, E$[d]=0$','GL-SigRep, E$[d]=10$','Interpreter','latex','Location','southeast')
grid on;
STANDARDIZE_FIGURE(fig)