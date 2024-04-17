clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
load("edge2.mat");
%% ER Curve Fit
figure;
subplot(1,3,1); hold on;
edge_spread_ER = edge_spread_ER/20^2;
plot(er_param,edge_spread_ER)
f = fit(er_param',edge_spread_ER,'poly3')
plot(er_param,f.p1*er_param.^3+ f.p2*er_param.^2 + f.p3*er_param + f.p4);
% plot(er_param,graphs.get_edge_spread('er',er_param,20));

%% BA Curve Fit
subplot(1,3,2); hold on;
edge_spread_BA = edge_spread_BA/20;
plot(ba_param,edge_spread_BA)
f = fit(ba_param',edge_spread_BA,'poly3')
plot(ba_param,f.p1*ba_param.^3+ f.p2*ba_param.^2 + f.p3*ba_param + f.p4);
% plot(ba_param,graphs.get_edge_spread('pa',ba_param,10));

%% RND Curve Fit
subplot(1,3,3); hold on;
edge_spread_RND = edge_spread_RND/20^2;
plot(rnd_param,edge_spread_RND)
f = fit(rnd_param',edge_spread_RND,'gauss1')
plot(rnd_param,f.a1*exp(-((rnd_param-f.b1)/f.c1).^2));
% plot(rnd_param,graphs.get_edge_spread('gaussian',rnd_param,20));
%%
load("edge2_n10.mat");
%% ER Curve Fit
subplot(1,3,1); hold on;
edge_spread_ER = edge_spread_ER/10^2;
plot(er_param,edge_spread_ER)
f = fit(er_param',edge_spread_ER,'poly3')
% plot(er_param,f.p1*er_param.^3+ f.p2*er_param.^2 + f.p3*er_param + f.p4);
% plot(er_param,graphs.get_edge_spread('er',er_param,20));

%% BA Curve Fit
subplot(1,3,2); hold on;
edge_spread_BA = edge_spread_BA/10;
plot(ba_param,edge_spread_BA)
f = fit(ba_param',edge_spread_BA,'poly3')
% plot(ba_param,f.p1*ba_param.^3+ f.p2*ba_param.^2 + f.p3*ba_param + f.p4);
% plot(ba_param,graphs.get_edge_spread('pa',ba_param,10));

%% RND Curve Fit
subplot(1,3,3); hold on;
edge_spread_RND = edge_spread_RND/10^2;
plot(rnd_param,edge_spread_RND)
f = fit(rnd_param',edge_spread_RND,'gauss1')
% plot(rnd_param,f.a1*exp(-((rnd_param-f.b1)/f.c1).^2));
% plot(rnd_param,graphs.get_edge_spread('gaussian',rnd_param,20));