clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('testmat.mat');
AR_params.N = 20;
AR_params.P = 20;
AR_params.gamma = 500;%2.6;
AR_params.threshold = .35;%.085;

N = 20;
A_e = zeros(N);
for i = 1:length(ce_i)
    A_e(ce_i(i)+1, ce_j(i)+1) = 1;
end
A_e = A_e+A_e';
L_0 = graphs.to_laplacian(A_e);

graphs.plot(A_e+A_e','A');
figure
graphs.plot(A_e,'A')

figure;
signals.plot(out(:,1400:1600)')
%%
gammas = logspace(-8,8,40);
thresholds = logspace(-4,-1,25);
t_max = length(thresholds);
trials = 20;
F = zeros(length(gammas),length(thresholds));
figure
for k = 1:length(gammas)
    toc
    tic
    k/length(gammas)*100
    AR_params.gamma = gammas(k);
    f = zeros(trials,1);
    L = GL.AR(out,AR_params);
    graphs.vcompare(L_0,L,'L_e','L')
    for t = 1:t_max 
        L_tmp = GL.threshold(L,thresholds(t));
        [~,~,f(t,i),~,~] = graphs.performance(L_0,L_tmp);
    end
    F(k,:) = mean(f,2);
end
%%
figure;
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
s = surf(X_axis,Y_axis,F');
