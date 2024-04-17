clc;clear;close all hidden;
figure;hold on;
load('large_performance_data_results_diff_n_er.mat')
F_list2 = F_list;
N_list2 = N_list;
load('large_performance_data_results_diff_n_er_100.mat')
F_list(end+1:end+length(F_list2)) = F_list2;
N_list(end+1:end+length(N_list2)) = N_list2;
plot(N_list,F_list,'LineWidth',2)

load('large_maximum_data_results_count_er.mat')
F_list2 = F_list;
N_list2 = N_list;
load('large_maximum_data_results_count_er_100.mat')
F_list(end+1:end+length(F_list2)) = F_list2;
N_list(end+1:end+length(N_list2)) = N_list2;
plot(N_list,F_list,'--','LineWidth',2);

load('large_maximum_data_results_count_er_100_nodelay.mat')
plot(N_list,F_list,'--','LineWidth',2);


load('large_maximum_data_results_count_er_100_newP.mat')
plot(N_list,F_list,'--');

load('large_maximum_data_results_count_dong_100_newP.mat')
plot(N_list,F_list,'LineWidth',2);

load('large_maximum_data_results_count_dong_100_newP_nodelay.mat')
N_list2 = N_list;
F_list2 = F_list;
load('large_maximum_data_results_count_dong_100_newP_nodelay2.mat')
N_list = [N_list N_list2];
F_list = [F_list;F_list2];
plot(N_list,F_list,'LineWidth',2);

load('large_maximum_data_results_count_er_100_undirected.mat')
plot(N_list,F_list,'LineWidth',2);

load('large_maximum_data_results_count_er_100_increasedSNR.mat')
plot(N_list,F_list,'LineWidth',2);
% 
load('large_maximum_data_results_count_er_100_increasedSNR2.mat')
plot(N_list,F_list,'LineWidth',2);


load('large_maximum_data_results_count_er_100_increasedSNR3.mat')
plot(N_list,F_list,'LineWidth',2);


load('large_maximum_data_results_count_er_100_increasedSNR4.mat')
plot(N_list,F_list,'LineWidth',2);

xlabel('N')
ylabel('Fscore')
title(signal_params.graph_type)
ylim([0 1])
grid on;

legend('AR algorithm','AR bruteforce','AR bruteforce w/ no delay','AR bruteforce w/ high AR','dong bruteforce', 'dong bruteforce w/ no delay','AR bruteforce /w undirected','AR bruteforce w/ increased SNR','AR bruteforce w/ increased SNR200','AR bruteforce w/ increased SNR200 (AR40)','AR bruteforce w/ increased SNR20000 (AR=200)')