clc;clear;close all hidden;
% load('large_maximum_data_results_n20_3.mat');
% load('large_maximum_data_results_n15_n10_3.mat');
% load('large_maximum_data_results_n40.mat');
load('large_new_maximum_data_results_n40_n20_n10.mat');

addpath(genpath('../../GL_classes/'));
gammas = logspace(-4,2,60);
thresholds = 0:.0025:2;
%%
figure;
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds(4:end)*signal_params.N);
s = surf(X_axis,Y_axis,F_ER_20(:,4:end,10)');
s.EdgeColor = "none";
xlabel("log$$_{10}(\gamma)$$",'Interpreter','latex','FontSize',14)
ylabel("$$\tilde{\delta}$$",'Interpreter','latex','FontSize',14)
zlabel("F-score",'Interpreter','latex','FontSize',14)
view(2)
%%
figure;hold on;
fig.fig = gcf;
%%
max_gammas = zeros(size(S_ER_20));
for i = 1:length(S_ER_20)
    [~, ind] = max(F_ER_20(:,4:end,i),[],"all");
    [row, col] = ind2sub(size(F_ER_20(:,:,i)),ind);
    max_gammas(i) = log10(gammas(row));
end
% fig.p1 = scatter(S_ER_20,(max_gammas),'filled','SizeData',100)
%%
max_gammas2 = zeros(size(S_ER_40));
for i = 1:length(S_ER_40)
    [~, ind] = max(F_ER_40(:,4:end,i),[],"all");
    [row, col] = ind2sub(size(F_ER_40(:,:,i)),ind);
    max_gammas2(i) = log10(gammas(row));
end
% plot(S_ER_40,(max_gammas2))
%%
max_gammas2 = zeros(size(S_ER_10));
for i = 1:length(S_ER_10)
    [~, ind] = max(F_ER_10(:,4:end,i),[],"all");
    [row, col] = ind2sub(size(F_ER_10(:,:,i)),ind);
    max_gammas2(i) = log10(gammas(row));
end
%%
fig.p2 = scatter([S_ER_10;S_ER_20],[max_gammas2;max_gammas],'filled','SizeData',100,'MarkerFaceColor','red')
%%
S = [S_ER_20;S_ER_10];
g = [max_gammas;max_gammas2];
%%
% a = -67.8256;
% b = 0.0413;
% c = 69.3318;
% plot(S_ER_20,c+a*exp(b.*S_ER_20.^2))

a = -0.2724;
b = 2.6540;
plot(S_ER_20,a*exp(b.*S_ER_20),'LineWidth',3,'color','blue');
%%
legend('20','40','10','alg');
legend('Empirical Data','Eq. 3.28','Interpreter','latex');
ylabel('Optimal log$_{10}(\gamma)$','Interpreter','latex')
xlabel('Graph Density, $\rho$','Interpreter','latex')
STANDARDIZE_FIGURE(fig)
grid on