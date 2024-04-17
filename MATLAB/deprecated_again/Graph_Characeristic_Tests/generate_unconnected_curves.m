% clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
load("unconnected_vertices.mat");
%% ER Curve Fit
figure;
subplot(1,3,1); hold on;
plot(er_param,unconnected_ER/signal_params.N)
% plot(er_param,graphs.get_unconnected_vertices('er',er_param));

%% BA Curve Fit
subplot(1,3,2); hold on;
plot(ba_param,unconnected_BA/signal_params.N)
% plot(er_param,graphs.get_unconnected_vertices('er',er_param));

%% RND Curve Fit
subplot(1,3,3); hold on;
plot(rnd_param,unconnected_RND/signal_params.N)
% plot(er_param,graphs.get_unconnected_vertices('er',er_param));