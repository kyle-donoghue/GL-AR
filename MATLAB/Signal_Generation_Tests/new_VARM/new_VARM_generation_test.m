clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
signal_params.N = 4;
signal_params.raw = 1;

signal_params.SNR = 2e3;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;

signal_params.intervals = 1;
signal_params.interval_length = 128;
signal_params.M = signal_params.intervals*signal_params.interval_length;

signal_params.order = 1;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .1;
signal_params.pole_mean = .7;
signal_params.delay = 0;

%%
x = zeros(signal_params.N,signal_params.M);
x(1,6:17) = sin(.25*(1:12));
figure(1);hold;
plot(x')

%%
A_0 = zeros(signal_params.N);
A_0(2,1) = 1;
D = graphs.create_VARM_tensor(signal_params,A_0);
a = permute(D(2,1,:),[3 1 2]);
y2 = filter(a(1),a(2:end),x(1,:));
plot(y2','--')

y = signals.customFilter(signal_params,x,D);
plot(y')
