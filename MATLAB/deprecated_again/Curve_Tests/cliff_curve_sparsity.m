%%
% clc;clear;close all hidden;
%%
load('gamma_threshold_sweep_with_sparsity_big.mat');
%%
F_t = .001;
%%
fits = zeros(length(sparsities),2);
for s = 1:length(sparsities)
    F = sparse_F(:,:,s);
    F_bin = double(F>F_t);
    y = zeros(length(gammas),1);
    for i = 1:length(gammas)
        indices = find(F_bin(i,:)==1);
        y(i) = thresholds(max(indices));
    end
    y = y(y<thresholds(end));
    y = y*signal_params.N-1;
    % y = y-1/signal_params.N;
    fit_data = fit(log10(gammas(1:length(y)))',y,'exp1');
    fits(s,:) = [fit_data.a fit_data.b];
end
fits(4,2) = 1.8;
%%
figure
subplot(1,2,1)
plot(sparsities,fits(:,1));
xlabel('Sparsity')
ylabel('a that produces correct CliffCurve')
title('a*e\^(bx)')
subplot(1,2,2)
plot(sparsities,fits(:,2));
xlabel('Sparsity')
ylabel('b that produces correct CliffCurve')
title('a*e\^(bx)')
%%
fit_data1 = fit(sparsities',fits(:,1)-fits(end,1),'exp1')
a_hat = fit_data1.a*exp(fit_data1.b*sparsities)+fits(end,1);
subplot(1,2,1);hold on;
plot(sparsities,a_hat)
%%
options = fitoptions('gauss1');
options.Lower(2) = 0;
options.Upper = [Inf 0 Inf];
offset = 1.45;
fit_data2 = fit(sparsities',fits(:,2)-offset,'gauss1',options)
b_hat = fit_data2.a1*exp(-((sparsities-fit_data2.b1)/fit_data2.c1).^2)+offset;
subplot(1,2,2);hold on;
plot(sparsities,b_hat)

