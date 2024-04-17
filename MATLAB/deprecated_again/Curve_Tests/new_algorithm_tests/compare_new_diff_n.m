clc;clear;close all hidden;
legend_list = {};
figure;hold on;
load('large_new_performance_data_results_diff_n_er.mat')
plot(N_list,F_list,'LineWidth',2)

load('large_maximum_data_results_count_er.mat')
F_list2 = F_list;
N_list2 = N_list;
load('large_maximum_data_results_count_er_100.mat')
F_list(end+1:end+length(F_list2)) = F_list2;
N_list(end+1:end+length(N_list2)) = N_list2;
plot(N_list,F_list,'--','LineWidth',2);

load("large_new_performance_data_results_diff_n_er_SNR200_AR100.mat")
plot(N_list,F_list,'LineWidth',2)

load('large_new_performance_data_results_diff_n_er_z.mat')
plot(N_list,F_list,'LineWidth',2)

load("large_new_performance_data_results_diff_n_er_z_SNR200_AR100.mat")
plot(N_list,F_list,'LineWidth',2)

load("large_new_performance_data_results_diff_n_er_z_SNR200_AR100_quadprog.mat")
plot(N_list,F_list,'LineWidth',2)

xlabel('N')
ylabel('Fscore')
title(signal_params.graph_type)
ylim([0 1])
grid on;

legend('AR algorithm','AR bruteforce','AR algorithm w/ SNR200,AR100','AR new algorithm','AR new algorithm w/ SNR200,AR100','AR new algorithm w/ SNR200,AR100 quadprog')