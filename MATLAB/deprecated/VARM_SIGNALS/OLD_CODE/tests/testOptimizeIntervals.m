clc;clear;close all hidden;
addpath('../creating_signals/');
addpath('../optimization/');
addpath('../accuracy/');
n = 12;
N = 4e3;
delay = 50;
interval = 10;
Fs = 250;
freqs = generate_frequencies(n,5,100,.1);
noise = .5;
P = 30;
[coeffs, y, x, A] = sineER(n,P,.4,delay,interval,N,Fs,freqs,noise);
alpha = 1;
beta = .15;

for i = 1:interval
    L_interval(:,:,i) = opt_L(coeffs(:,:,i)',alpha,beta);
end
L = mean(L_interval,3);

A>0
A2 = laplacian_to_adjacency(L)
figure
subplot(1,2,1)
imagesc((A>0))
subplot(1,2,2)
imagesc((A2>.09))
A = A>0;
threshold_array = .01:.01:.5;
for i = 1:length(threshold_array)
    threshold = threshold_array(i);
    A3 = A2>threshold;
    f(i) = calculate_F_score(A>0,A3);
end
for i = 1:interval
    A_interval(:,:,i) = laplacian_to_adjacency(L_interval(:,:,i));
end
for i = 1:length(threshold_array)
    threshold = threshold_array(i);
    A3 = mean(A_interval>threshold,3);
    f2(i) = calculate_F_score(A>0,A3);
end
figure;
subplot(1,2,1)
plot(threshold_array,f)
subplot(1,2,2)
plot(threshold_array,f2)