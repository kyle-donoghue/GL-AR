clc;clear;close all hidden;
load("large_new_performance_data_results_n20_n10.mat")
load("large_new_performance_data_results_n10.mat")
load("large_new_performance_data_results_n10_2.mat")
load("large_new_performance_data_results_n20_2.mat")

load("large_new_performance_data_results_n20_4.mat")
load("large_new_performance_data_results_n20_n10_3.mat")
load("large_new_performance_data_results_n10_4.mat")
load("large_new_performance_data_results_n20_n10_4.mat")
load("large_new_performance_data_results_n20_n10_5.mat")

%% N = 20
plot_perf(1,1,S_ER_20,F_ER_20);
plot_perf(1,2,S_BA_20,F_BA_20);
plot_perf(1,3,S_RND_20,F_RND_20);
%% N = 10
plot_perf(2,1,S_ER_10,F_ER_10);
plot_perf(2,2,S_BA_10,F_BA_10);
plot_perf(2,3,S_RND_10,F_RND_10);
%%
load("large_new_maximum_data_results_n40_n20_n10.mat");
%% N = 20
plot_max(1,1,S_ER_20,F_ER_20);
plot_max(1,2,(1:10)*2/20,F_BA_20(:,:,1:10));
plot_max(1,3,S_RND_20,F_RND_20);
%% N = 10
plot_max(2,1,S_ER_10,F_ER_10);
plot_max(2,2,(1:5)*2/10,F_BA_10(:,:,1:5));
plot_max(2,3,S_RND_10,F_RND_10);

%%
% figure;hold on;
% plot(S_RND_10,T_RND_10(:,1))
%  plot((1:5)*2/10,T_BA_10(:,1))
% plot(S_ER_10,T_ER_10(:,1))
% 
% figure;hold on;
% plot(S_RND_20,T_RND_20(:,1))
%  plot(S_BA_20,T_BA_20(:,1))
% plot(S_ER_20,T_ER_20(:,1))
% 
% scalar = 1;
% figure;
% subplot(1,3,1);hold on;
% plot(S_ER_10,T_ER_10(:,1),'--')
% plot(S_ER_10,T_ER_10(:,2)*scalar)
% subplot(1,3,2);hold on;
% plot((1:5)*2/10,T_BA_10(:,1),'--')
% plot((1:5)*2/10,T_BA_10(:,2)*scalar)
% subplot(1,3,3);hold on;
% plot(S_RND_10,T_RND_10(:,1),'--')
% plot(S_RND_10,T_RND_10(:,2)*scalar)

function plot_perf(i,j,s,f)
    figure(i);
    subplot(1,3,j);hold on;
    plot(s,f*100,'LineWidth',2);
end
function plot_max(i,j,s,f)
    figure(i);
    subplot(1,3,j);hold on;
    plot(s,permute(max(f,[],[1 2]),[3 1 2])*100,'--','LineWidth',2);
    xlabel('Graph Density, $$\rho$$','Interpreter','latex','FontSize',14)
    ylabel('F-score','Interpreter','latex','FontSize',14)
    % title('ER','Interpreter','latex','FontSize',14)
    legend('Algorithm Performance','Theoretical Performance')
    ylim([0 109])
end