%%
clc;clear;close all hidden;
%%
load('gamma_threshold_plot_for_cliff_fit.mat');
%%
F_t = .02;
%%
figure;
F_bin = double(F>F_t);
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
s = surf(X_axis,Y_axis,F_bin');
s.EdgeColor = "none";
%%
y = zeros(length(gammas),1);
for i = 1:length(gammas)
    indices = find(F_bin(i,:)==1);
    y(i) = thresholds(max(indices));
end
%%
figure;
s = surf(X_axis,Y_axis,F');
s.EdgeColor = "none";
hold on;
plot3(log10(gammas),y,.7*ones(length(gammas),1),LineWidth=4,color='red');
%%
y = y(y<thresholds(end));
fit_data = fit(log10(gammas(1:length(y)))',y-y(1),'exp1');
y_hat = fit_data.a*exp(fit_data.b*log10(gammas(1:length(y))));
plot3(log10(gammas(1:length(y))),y_hat+y(1),.7*ones(length(y),1),LineWidth=4,color='green');
%%
% fit_data.b = fit_data.b/1.25;
% fit_data.a = fit_data.a/2;
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1));
plot3(log10(gammas),x_hat+y(1),.7*ones(length(gammas),1),LineWidth=4,color='magenta');
%%
load('gamma_threshold_plot_for_cliff_fit2.mat');
%%
F_bin = double(F>F_t);
y = zeros(length(gammas),1);
for i = 1:length(gammas)
    indices = find(F_bin(i,:)==1);
    y(i) = thresholds(max(indices));
end

%%
figure;
s = surf(X_axis,Y_axis,F');
s.EdgeColor = "none";
hold on;
%%
y = y(y<thresholds(end));
y_hat = fit_data.a*exp(fit_data.b*log10(gammas(1:length(y))));
plot3(log10(gammas(1:length(y))),y_hat+y(1),.7*ones(length(y),1),LineWidth=4,color='green');
%%
% x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1));
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-.75));
x_hat = x_hat(x_hat<thresholds(end));
plot3(log10(gammas(1:length(x_hat))),x_hat+.95*y(1),.7*ones(length(x_hat),1),LineWidth=4,color='magenta');
%%
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1));
x_hat = x_hat(x_hat<thresholds(end));
plot3(log10(gammas(1:length(x_hat))),x_hat+.95*y(1),.7*ones(length(x_hat),1),LineWidth=4,color='magenta');
%%
load('gamma_threshold_plot_for_cliff_fit3.mat');
%%
F_t = .02;
F_bin = double(F>F_t);
y = zeros(length(gammas),1);
for i = 1:length(gammas)
    indices = find(F_bin(i,:)==1);
    y(i) = thresholds(max(indices));
end
%%
figure;
s = surf(X_axis,Y_axis,F_bin');
s.EdgeColor = "none";
%%
figure;
s = surf(X_axis,Y_axis,F');
s.EdgeColor = "none";
hold on;
%%
y = y(y<thresholds(end));
y_hat = fit_data.a*exp(fit_data.b*log10(gammas(1:length(y))));
plot3(log10(gammas(1:length(y))),y_hat+y(1),.7*ones(length(y),1),LineWidth=4,color='green');
fit_data = fit(log10(gammas(1:length(y)))',y-y(1),'exp1');
y_hat = fit_data.a*exp(fit_data.b*log10(gammas(1:length(y))));
plot3(log10(gammas(1:length(y))),y_hat+y(1),.7*ones(length(y),1),LineWidth=4,color='black');
%%
% x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1));
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-.75));
x_hat = x_hat(x_hat<thresholds(end));
plot3(log10(gammas(1:length(x_hat))),x_hat+.95*y(1),.7*ones(length(x_hat),1),LineWidth=4,color='magenta');
%%
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1));
x_hat = x_hat(x_hat<thresholds(end));
plot3(log10(gammas(1:length(x_hat))),x_hat+.95*y(1),.7*ones(length(x_hat),1),LineWidth=4,color='magenta');
%%
load('gamma_threshold_plot_for_cliff_fit4.mat');
%%
F_t = .02;
F_bin = double(F>F_t);
y = zeros(length(gammas),1);
for i = 1:length(gammas)
    indices = find(F_bin(i,:)==1);
    y(i) = thresholds(max(indices));
end
%%
figure;
s = surf(X_axis,Y_axis,F_bin');
s.EdgeColor = "none";
%%
figure;
s = surf(X_axis,Y_axis,F');
s.EdgeColor = "none";
hold on;
%%
y = y(y<thresholds(end));
y_hat = fit_data.a*exp(fit_data.b*log10(gammas(1:length(y))));
plot3(log10(gammas(1:length(y))),y_hat+y(1),.7*ones(length(y),1),LineWidth=4,color='green');
fit_data = fit(log10(gammas(1:length(y)))',y-y(1),'exp1');
y_hat = fit_data.a*exp(fit_data.b*log10(gammas(1:length(y))));
plot3(log10(gammas(1:length(y))),y_hat+y(1),.7*ones(length(y),1),LineWidth=4,color='black');
%%
% x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1));
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-.75));
x_hat = x_hat(x_hat<thresholds(end));
plot3(log10(gammas(1:length(x_hat))),x_hat+.95*y(1),.7*ones(length(x_hat),1),LineWidth=4,color='magenta');
%%
x_hat = fit_data.a*exp(fit_data.b*(log10(gammas)-1.5));
x_hat = x_hat(x_hat<thresholds(end));
plot3(log10(gammas(1:length(x_hat))),x_hat+.95*y(1),.7*ones(length(x_hat),1),LineWidth=4,color='magenta');
