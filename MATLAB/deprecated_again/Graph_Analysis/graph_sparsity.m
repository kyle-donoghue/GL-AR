clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
signal_params.N = 20;
signal_params.raw = 1;

%%
trials = 100;
sparsity = zeros(trials,3);
for i = 1:trials;
    [L_01,~] = graphs.create(signal_params,'er',.2);
    [L_02,~] = graphs.create(signal_params,'pa',1);
    [L_03,~] = graphs.create(signal_params,'gaussian',0.75,0.5);
    [~,~,~,~,n_e1] = graphs.performance(L_01,L_01);
    [~,~,~,~,n_e2] = graphs.performance(L_02,L_02);
    [~,~,~,~,n_e3] = graphs.performance(L_03,L_03);
    sparsity(i,1) = n_e1/(signal_params.N^2-signal_params.N);
    sparsity(i,2) = n_e2/(signal_params.N^2-signal_params.N);
    sparsity(i,3) = n_e3/(signal_params.N^2-signal_params.N);
end
sparsity = mean(sparsity)
