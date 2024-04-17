clc;clear;close all hidden;
addpath('../creating_signals/');
n = 4;
N = 10e3;
delay = 50;
interval = 1;
Fs = 200;
freqs = generate_frequencies(4,10,80)
noise = .2;
P = 10;
[coeffs, y, x, A] = sineER(n,P,.4,delay,interval,N,Fs,freqs,noise);
figure;
subplot(1,2,1);
plot(x);
subplot(1,2,2);
plot(y);
figure;
for i = 1:n
    subplot(2,4,i);
    plot(abs(fft(x(6000:7000,i))))
    subplot(2,4,i+4);
    plot(abs(fft(y(6000:7000,i))))
end
figure;
for i = 1:n
    subplot(2,4,i);
    plot(coeffs(:,i))
    subplot(2,4,i+4);
    plot(abs(freqz(1,coeffs(:,i))))
end
A
figure;
subplot(1,2,1)
periodogram(x(:,1))
subplot(1,2,2)
periodogram(y(:,1))