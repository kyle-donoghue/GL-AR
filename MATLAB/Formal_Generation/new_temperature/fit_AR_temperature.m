clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('temperature_dataset.mat')
%%
y_noisy = signals.z_score(detrended(:,1:2:end));
[fits,P] = GL.sweep_fits(y_noisy);
plot(P,fits)