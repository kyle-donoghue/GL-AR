clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
rng("default")
%%
signal_params.N = 8;
signal_params.raw = 1;

signal_params.intervals = signal_params.N;
signal_params.interval_length = 64;
signal_params.M = signal_params.intervals*signal_params.interval_length;

signal_params.order = 1;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .25;
signal_params.pole_mean = .5;

%%
[~,A] = graphs.create(signal_params,'er',.2);
A_0 = graphs.to_directed(A);

figure;
graphs.vcompare(A,A_0,'A','A_0');

%%
G = graphs.createGraphTensor(signal_params,A_0);
figure;
plot(abs(permute(G(2,3,:),[3 1 2])));
%%
x = randn(signal_params.N,signal_params.M);
x = zeros(signal_params.N,signal_params.M);
x(3,:) = sin(.2*(1:signal_params.M));
x(2,:) = sin(.08*(1:signal_params.M));
figure;plot(x')
X = signals.createTensorSpectrum(signal_params,x);

%%
Y = signals.filterTensorSpectrum(signal_params,X,G);
y = signals.inverseTensorSpectrum(signal_params,Y);
figure;plot(y')