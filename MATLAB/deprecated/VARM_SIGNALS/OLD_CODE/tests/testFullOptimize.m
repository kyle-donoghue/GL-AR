clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
n = 20;
N = 4e3;
delay = 50;
interval = 10;
Fs = 250;
% freqs = generate_frequencies(n,5,1000,.1)/10;
freqs = generate_frequencies2(n,5,100);
noise = .1;
P = 20;
[coeffs, y, x, A] = sineER(n,P,.4,delay,interval,N,Fs,freqs,noise);
%%
alpha = 1;
beta = .15;
tolerance = 1e-3;
solve_intervals = zeros(interval,1);
for i = 1:interval
    [L_interval(:,:,i),solve_intervals(i)] = full_optimization(coeffs(:,:,i)',alpha,beta,tolerance);
end
L = mean(L_interval,3);
for i = 1:interval
    % [L_interval2(:,:,i),solve_intervals(i)] = full_optimization(coeffs(2:end,:,i)',alpha,beta,tolerance);
    [L_interval2(:,:,i)] = opt_L(coeffs(:,:,i)',alpha,beta);
end
L2 = mean(L_interval2,3);
%%
A>0
A2 = laplacian_to_adjacency(L)
A4 = laplacian_to_adjacency(L2)

A = A>0;
threshold_array = .01:.01:.5;
for i = 1:length(threshold_array)
    threshold = threshold_array(i);
    A3 = A2>threshold;
    A5 = A4>threshold;
    error = (A-A3).^2;
    e(i) = mean(error(:));
    f(i) = calculate_F_score(A>0,A3);
    f2(i) = calculate_F_score(A>0,A5);
end
figure;
plot(threshold_array,e)
hold on;
plot(threshold_array,f,'--')
plot(threshold_array,f2,'--')
legend('e','f1','f2')
max(f)

thresholds1 = threshold_array(find(f==max(f)));
thresholds2 = threshold_array(find(f2==max(f2)));


figure
subplot(1,3,1)
imagesc((A>0))
subplot(1,3,2)
imagesc((A2>thresholds1(1)))
subplot(1,3,3)
imagesc((A4>thresholds2(1)))