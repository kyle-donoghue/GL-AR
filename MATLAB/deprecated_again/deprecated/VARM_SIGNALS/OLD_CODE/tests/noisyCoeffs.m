clc;clear;close all hidden;
addpath('../creating_signals/');
n = 1;
N = 400;
delay = 50;
interval = 1;
Fs = 250;
freqs = 20;
P = 100;
noise1 = .1;
noise2 = 5;
coeffs1 = zeros(P+1,n);
coeffs2 = coeffs1;
m = 25;
for i = 1:m
    [coeffs, y,x,A] = sineER(n,P,.4,delay,interval,N,Fs,freqs,noise1);
    coeffs1 = coeffs1+coeffs;
    [coeffs, y2, x2, A2] = sineER(n,P,.4,delay,interval,N,Fs,freqs,noise2);
    coeffs2 = coeffs2+coeffs;
end
coeffs1 = coeffs1/m;
coeffs2 = coeffs2/m;
figure
subplot(1,2,1);
plot(coeffs1)
subplot(1,2,2);
plot(coeffs2)
figure
freqz(1,coeffs1)
figure
freqz(1,coeffs2)