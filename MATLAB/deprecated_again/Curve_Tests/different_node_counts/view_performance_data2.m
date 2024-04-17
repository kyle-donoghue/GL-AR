clc;clear;close all hidden;
% load("large_performance_data_results_n20_4.mat")
% load("large_performance_data_results_n15_n10_4.mat")
% load("large_performance_data_results_n15_n10_5.mat")
% load("large_performance_data_results_n20_n15_n10_5.mat")
load("large_performance_data_results_n20_n15_n10_6.mat")
%% N = 10
figure(1);
subplot(1,3,1);hold on;
plot(S_ER_10,F_ER_10*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('ER','Interpreter','latex','FontSize',14)
subplot(1,3,2);hold on;
plot(S_BA_10,F_BA_10*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('BA','Interpreter','latex','FontSize',14)
subplot(1,3,3);hold on;
plot(S_RND_10,F_RND_10*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('RND','Interpreter','latex','FontSize',14)
%% N = 15
figure(2);
subplot(1,3,1);hold on;
plot(S_ER_15,F_ER_15*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('ER','Interpreter','latex','FontSize',14)
subplot(1,3,2);hold on;
plot(S_BA_15,F_BA_15*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('BA','Interpreter','latex','FontSize',14)
subplot(1,3,3);hold on;
plot(S_RND_15,F_RND_15*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('RND','Interpreter','latex','FontSize',14)
%% N = 20
figure(3);
subplot(1,3,1);hold on;
plot(S_ER_20,F_ER_20*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('ER','Interpreter','latex','FontSize',14)
subplot(1,3,2);hold on;
plot(S_BA_20,F_BA_20*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('BA','Interpreter','latex','FontSize',14)
subplot(1,3,3);hold on;
plot(S_RND_20,F_RND_20*100,'LineWidth',2);
xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
ylabel('F-score','Interpreter','latex','FontSize',14)
title('RND','Interpreter','latex','FontSize',14)
%%
load("large_maximum_data_results_n20_3.mat")
load("large_maximum_data_results_n15_n10_3.mat")
%% N = 10
figure(1);
subplot(1,3,1);hold on;
plot(S_ER_10,permute(max(F_ER_10,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
subplot(1,3,2);hold on;
plot(S_BA_10,permute(max(F_BA_10,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
subplot(1,3,3);hold on;
plot(S_RND_10,permute(max(F_RND_10,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
%% N = 15
figure(2);
subplot(1,3,1);hold on;
plot(S_ER_15,permute(max(F_ER_15,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
subplot(1,3,2);hold on;
plot(S_BA_15,permute(max(F_BA_15,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
subplot(1,3,3);hold on;
plot(S_RND_15,permute(max(F_RND_15,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
%% N = 20
figure(3);
subplot(1,3,1);hold on;
plot(S_ER_20,permute(max(F_ER_20,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
subplot(1,3,2);hold on;
plot(S_BA_20,permute(max(F_BA_20,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
subplot(1,3,3);hold on;
plot(S_RND_20,permute(max(F_RND_20,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
legend('Algorithm Performance','Theoretical Performance','Interpeter','latex')
ylim([0 109])
