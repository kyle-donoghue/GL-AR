%%
clc;clear;close all hidden;
%%
load('gamma_sparsity_relationship_v3.mat');
%%
figure;hold on;
plot(log10(gammas),gamma_curves')

%%
new_max_gamma = max_gamma;
% new_max_gamma(new_max_gamma < 10^(-.12)) = 10^(-.12);
new_max_gamma(new_max_gamma < 10^(-.3)) = 10^(-.4);
%%
figure;hold on;
plot(sparsities,log10(max_gamma))
plot(sparsities,log10(new_max_gamma))
%%
% plot(sparsities,.5./sparsities.^2)
% plot(sparsities,100*exp(-sparsities.^2/.0255))
%%

options = fitoptions('logistic4')
options.Lower = [0 -Inf -Inf -Inf];
options.Upper = [0 Inf Inf Inf];
% fit_data = fit(sparsities',log10(new_max_gamma)+.12,'logistic4',options)
fit_data = fit(sparsities',log10(new_max_gamma)+.4,'logistic4',options)
y_hat = fit_data.d+(0-fit_data.d)./(1+(sparsities./fit_data.c).^fit_data.b);
figure;hold on;
plot(sparsities,log10(new_max_gamma)+.4)
plot(sparsities,y_hat);

%%
figure;
for i = 1:length(sparsities)
    plot(log10(gammas),gamma_curves(i,:));hold on;
    xline(y_hat(i)-.4);hold off;
    title(sparsities(i))
    pause(1)
end