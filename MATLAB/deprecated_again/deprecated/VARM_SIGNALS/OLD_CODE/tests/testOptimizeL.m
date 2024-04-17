clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
n = 8;
N = 4e3;
delay = 50;
interval = 1;
Fs = 250;
freqs = generate_frequencies(n,5,100,.1);
noise = .2;
P = 20;
[coeffs, y, x, A] = sineER(n,P,.4,delay,interval,N,Fs,freqs,noise);
alpha = 1;
beta = .1;
L = opt_L(coeffs',alpha,beta);
A>0
A2 = laplacian_to_adjacency(L)
figure
subplot(1,2,1)
imagesc((A>0))
subplot(1,2,2)
imagesc((A2>.15))
A = A>0;
P = calculate_precision(A>0,A2>.15)
R = calculate_recall(A>0,A2>.15)
F = calculate_F_score(A>0,A2>.15)
[tp,tn,fp,fn] = calculate_result_true_false(A>0,A2>.15);
threshold_array = .01:.01:.5;
for i = 1:length(threshold_array)
    threshold = threshold_array(i);
    A3 = A2>threshold;
    error = (A-A3).^2;
    e(i) = mean(error(:));
    f(i) = calculate_F_score(A>0,A2>threshold);
end
figure;
plot(threshold_array,e);
hold on;
plot(threshold_array,f,'--')