clc;clear;close all hidden;
load("large_performance_data_results_n20_4.mat")
% load("large_performance_data_results_n40.mat")
% %% N = 10
% figure(1);
% subplot(1,3,1);hold on;
% plot(S_ER_10,F_ER_10*100);
% subplot(1,3,2);hold on;
% plot(S_BA_10,F_BA_10*100);
% subplot(1,3,3);hold on;
% plot(S_RND_10,F_RND_10*100);
%% N = 20
figure(2);
subplot(1,3,1);hold on;
plot(S_ER_20,F_ER_20*100);
subplot(1,3,2);hold on;
plot(S_BA_20,F_BA_20*100);
subplot(1,3,3);hold on;
plot(S_RND_20,F_RND_20*100);
% %% N = 40
% figure(3);
% subplot(1,3,1);hold on;
% plot(S_ER_40,F_ER_40*100);
% subplot(1,3,2);hold on;
% plot(S_BA_40,F_BA_40*100);
% subplot(1,3,3);hold on;
% plot(S_RND_40,F_RND_40*100);
%%
load("large_maximum_data_results_n20_3.mat")
% load("large_maximum_data_results_n15_n10_3.mat")

% load("large_maximum_data_results_n20_2.mat")
% load("large_maximum_data_results_n40.mat")
% %% N = 10
% figure(1);
% subplot(1,3,1);hold on;
% plot(S_ER_10,permute(max(F_ER_10,[],[1 2]),[3 1 2])*100,'--');
% ylim([0 100])
% subplot(1,3,2);hold on;
% plot(S_ER_10,permute(max(F_ER_10,[],[1 2]),[3 1 2])*100,'--');
% ylim([0 100])
% subplot(1,3,3);hold on;
% plot(S_ER_10,permute(max(F_ER_10,[],[1 2]),[3 1 2])*100,'--');
% ylim([0 100])
%% N = 20
figure(2);
subplot(1,3,1);hold on;
plot(S_ER_20,permute(max(F_ER_20,[],[1 2]),[3 1 2])*100,'--');
ylim([0 100])
subplot(1,3,2);hold on;
plot(S_BA_20,permute(max(F_BA_20,[],[1 2]),[3 1 2])*100,'--');
ylim([0 100])
subplot(1,3,3);hold on;
plot(S_RND_20,permute(max(F_RND_20,[],[1 2]),[3 1 2])*100,'--');
ylim([0 100])
% %% N = 40
% figure(3);
% subplot(1,3,1);hold on;
% plot(S_ER_40,permute(max(F_ER_40,[],[1 2]),[3 1 2])*100,'--');
% ylim([0 100])
% subplot(1,3,2);hold on;
% plot(S_BA_40,permute(max(F_BA_40,[],[1 2]),[3 1 2])*100,'--');
% ylim([0 100])
% subplot(1,3,3);hold on;
% plot(S_RND_40,permute(max(F_RND_40,[],[1 2]),[3 1 2])*100,'--');
% ylim([0 100])
