clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
freq = 20;
Fs = 100;
x = 10*sin(2*pi*freq*(1/Fs:1/Fs:1000/Fs));
x = x+0*randn(size(x));
plot(x)
%%
% A = [1 -2.7607 3.8106 -2.6535 0.9238];
% x = randn(1000,1);
% x = filter(1,A,x)';
%%
% [fits,P] = GL.sweep_fits(x,50);
% figure
% plot(P,fits)
%%
P = 200;
b = signals.approxAR(x,P);
[H_P,w_P] = periodogram(x,[],1024);
[H_AR,w_AR] = freqz(1,[1 b],1024);
H = (abs(H_AR)-1)/max(abs(H_AR))*max(sqrt(H_P));
figure;hold on;
fig.fig = gcf;
% plot(w_AR,abs(H_AR))
fig.p1 = plot(w_P,10*log10(H_P),'--','LineWidth',1);
fig.p2 = plot(w_AR,10*log10(H.^2),'linewidth',1);
xlabel('Normalized Radians','Interpreter','latex')
ylabel('Power Spectrum (dB)','Interpreter','latex')
legend('Periodogram  ','AR Spectrum  ','Interpreter','latex')
STANDARDIZE_FIGURE(fig)
grid on
%%
% y = sin([1:300]') + 0.5*randn(300,1);
% y = iddata(y);
% sys_fb = ar(y,4,'fb');
% figure
% spectrum(sys_fb)
