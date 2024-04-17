%%
close all hidden;
load('weight_distribution_500_trials2.mat');
%%
figure;
for i = length(gammas):-1:1
    mean(big_weights2(i,:))
    histogram(big_weights2(i,:),n_hist);
    th = .1;
    xline(th);
    %xlim([0 1]); % uncomment this to fix x axis
    title(sprintf("Gamma=%d, F=%d",gammas(i),f(i,find(thresholds==th))));
    % pause(.5)
end
%%
n_hist = 50;
figure(2);
figure(3);
pause(2);
for i = 1:length(gammas)
    mean(big_weights2(i,:))
    figure(2);
    histogram(big_weights2(i,:),n_hist);
    [~,ind] = max(f(i,:));
    xline(thresholds(ind));
    % xlim([0 .5]); % uncomment this to fix x axis
    title(sprintf("Gamma=%d, F=%d",i,f(i,ind)));
    figure(3);
    hold off;surf(f');hold on; scatter3(i,ind,f(i,ind),100,'red','filled')
    pause(1.5)
end
%%
n_hist = 25;
figure
subplot(1,3,1)
histogram(big_weights2(20,:),n_hist);
xline(thresholds(12));
title(sprintf("Gamma=%.1e, F=%.1d",gammas(20),f(20,12)*100));
subplot(1,3,2)
histogram(big_weights2(11,:),n_hist);
xline(thresholds(12));
title(sprintf("Gamma=%.1e, F=%.1d",gammas(15),f(11,12)*100));
subplot(1,3,3)
histogram(big_weights2(3,:),n_hist);
xline(thresholds(12));
title(sprintf("Gamma=%.1e, F=%.1d",gammas(3),f(3,12)*100));

%% README
% I think that once the distribution turns gaussian, but is wide enough for
% the threshold to still be on the distribution -> where F_max occurs.

% Maybe F_max occurs when the threshold is at the point where it equals
% sparsity of the normal dist, like related to Z-score.

% For example, it would be 2*[CDF_gauss(threshold)-.5]=sparsity, or
% something along those lines