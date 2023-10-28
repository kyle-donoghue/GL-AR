clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
signal_params.N = 4;
signal_params.raw = 1;

signal_params.SNR = 2e0;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;

signal_params.intervals = signal_params.N;
signal_params.interval_length = 512;
signal_params.M = signal_params.intervals*signal_params.interval_length;

signal_params.order = 1;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .25;
signal_params.pole_mean = .5;

signal_params.minSeparation = 10;
signal_params.maxSeparation = 100;
%%
[L_0,A] = graphs.create(signal_params,'er',.4);
A_0 = graphs.to_directed(A)

G = graphs.createGraphTensor(signal_params,A_0);

x = randn(signal_params.N,signal_params.M);
x = signals.create_raw_sine(signal_params);
x = signals.createFullSinePulse(signal_params);


X = signals.createTensorSpectrum(signal_params,x);

Y = signals.filterTensorSpectrum(signal_params,X,G);

y = signals.inverseTensorSpectrum(signal_params,Y);
%%
figure;
subplot(1,2,1);
signals.plot(x');
subplot(1,2,2);
signals.plot(y');

%%
graphs.view_tensor_spectrum(G)

