clc;clear;close all hidden;
addpath('../creating_signals/');
n = 4;
N = 4e3;
delay = 200;
interval = 1;
innov = 1e3;
Fs = 1e3;
freqs = [2;3;5;7];
noise = 0;
P = 20;
[coeffs, y, x, A] = sineBA(n,P,delay,interval,innov,N,Fs,freqs,noise);
figure;
subplot(1,2,1);
plot(x);
subplot(1,2,2);
plot(y);