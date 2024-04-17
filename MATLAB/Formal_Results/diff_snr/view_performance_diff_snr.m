clc;clear;close all hidden;
addpath('../plot_functions');
%%
figure;
fig.fig = gcf;
%%
load('../../Formal_Data/diff_snr/performance_data_diff_snr_sine_FINAL_rect.mat')
fig.p1 = plot_perf(1,1,10*log10(snr_list),F_ER_20(:,3),3);
fig.p2 = plot_perf(1,2,10*log10(snr_list),F_BA_20(:,3),3);
fig.p3 = plot_perf(1,3,10*log10(snr_list),F_RND_20(:,3),3);
%%
load('../../Formal_Data/diff_snr/maximum_data_diff_snr_dong_sine_FINAL_nodelay_rect.mat');
fig.p4 = plot_max(1,1,10*log10(snr_list),F_ER_surfs,3);
fig.p5 = plot_max(1,2,10*log10(snr_list),F_BA_surfs,3);
fig.p6 = plot_max(1,3,10*log10(snr_list),F_RND_surfs,3);
%%
load('../../Formal_Data/diff_snr/maximum_data_diff_snr_dong_sine_FINAL_10delay_rect.mat');
fig.p7 = plot_max(1,1,10*log10(snr_list),F_ER_surfs,3);
fig.p8 = plot_max(1,2,10*log10(snr_list),F_BA_surfs,3);
fig.p9 = plot_max(1,3,10*log10(snr_list),F_RND_surfs,3);
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