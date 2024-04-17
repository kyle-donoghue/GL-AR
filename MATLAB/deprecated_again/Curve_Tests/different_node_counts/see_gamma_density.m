clc;clear;close all hidden;
load('large_maximum_data_results_n20_3.mat');
addpath(genpath('../../GL_classes/'));
gammas = logspace(-2,4,60);
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
%%
max_gammas = zeros(size(S_ER_20));
for i = 1:length(S_ER_20)
    [~, ind] = max(F_ER_20(:,4:end,i),[],"all");
    [row, col] = ind2sub(size(F_ER_20(:,:,i)),ind);
    max_gammas(i) = log10(gammas(row));
end
plot(S_ER_20,(max_gammas))
%%
max_gammas2 = zeros(size(S_RND_20));
for i = 1:length(S_RND_20)
    [~, ind] = max(F_RND_20(:,4:end,i),[],"all");
    [row, col] = ind2sub(size(F_RND_20(:,:,i)),ind);
    max_gammas2(i) = log10(gammas(row));
end
plot(S_RND_20,(max_gammas2))
%%
S = [S_ER_20;S_RND_20];
g = [max_gammas;max_gammas2];
%%
a = -67.8256;
b = 0.0413;
c = 69.3318;
plot(S_RND_20,c+a*exp(b.*S_RND_20.^2))
