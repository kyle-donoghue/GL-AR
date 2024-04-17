% clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
load("variance_vertices.mat");
%% ER Curve Fit
figure;
subplot(1,3,1); hold on;
plot(er_param,std_ER)
plot(er_param,graphs.get_vertice_std('er',er_param,signal_params.N));
o = fitoptions('poly3');
f = fit(er_param',std_ER,'poly3',o)
plot(er_param,f.p1*er_param.^3+ f.p2*er_param.^2 + f.p3*er_param + f.p4);
%% BA Curve Fit
subplot(1,3,2); hold on;
plot(ba_param,std_BA)
plot(ba_param,graphs.get_vertice_std('pa',ba_param,signal_params.N));
o = fitoptions('poly3');
f = fit(ba_param',std_BA,'poly3',o)
plot(ba_param,f.p1*ba_param.^3+ f.p2*ba_param.^2 + f.p3*ba_param + f.p4);
%% RND Curve Fit
subplot(1,3,3); hold on;
plot(rnd_param,std_RND)
plot(rnd_param,graphs.get_vertice_std('gaussian',rnd_param,signal_params.N));
o = fitoptions('poly3');
f = fit(rnd_param',std_RND,'poly3',o)
plot(rnd_param,f.p1*rnd_param.^3+ f.p2*rnd_param.^2 + f.p3*rnd_param + f.p4);
%%
load("variance_vertices_n10.mat");
ratio = 20/10;
%% ER Curve Fit
subplot(1,3,1); hold on;
plot(er_param,std_ER*sqrt(ratio))
% plot(er_param,graphs.get_vertice_std('er',er_param,signal_params.N));
o = fitoptions('poly3');
f = fit(er_param',std_ER,'poly3',o)
% plot(er_param,f.p1*er_param.^3+ f.p2*er_param.^2 + f.p3*er_param + f.p4);
%% BA Curve Fit
subplot(1,3,2); hold on;
plot(ba_param,std_BA*sqrt(ratio))
% plot(ba_param,graphs.get_vertice_std('pa',ba_param,signal_params.N));
o = fitoptions('poly3');
f = fit(ba_param',std_BA,'poly3',o)
% plot(ba_param,f.p1*ba_param.^3+ f.p2*ba_param.^2 + f.p3*ba_param + f.p4);
%% RND Curve Fit
subplot(1,3,3); hold on;
plot(rnd_param,std_RND*max(1.8,sqrt(ratio)))
% plot(rnd_param,graphs.get_vertice_std('gaussian',rnd_param,signal_params.N));
o = fitoptions('poly3');
f = fit(rnd_param',std_RND,'poly3',o)
% plot(rnd_param,f.p1*rnd_param.^3+ f.p2*rnd_param.^2 + f.p3*rnd_param + f.p4);
%%
load("variance_vertices_n15.mat");
ratio = 20/15;
%% ER Curve Fit
subplot(1,3,1); hold on;
plot(er_param,std_ER*sqrt(ratio))
% plot(er_param,graphs.get_vertice_std('er',er_param,signal_params.N));
o = fitoptions('poly3');
f = fit(er_param',std_ER,'poly3',o)
% plot(er_param,f.p1*er_param.^3+ f.p2*er_param.^2 + f.p3*er_param + f.p4);
%% BA Curve Fit
subplot(1,3,2); hold on;
plot(ba_param,std_BA*sqrt(ratio))
% plot(ba_param,graphs.get_vertice_std('pa',ba_param,signal_params.N));
o = fitoptions('poly3');
f = fit(ba_param',std_BA,'poly3',o)
% plot(ba_param,f.p1*ba_param.^3+ f.p2*ba_param.^2 + f.p3*ba_param + f.p4);
%% RND Curve Fit
subplot(1,3,3); hold on;
plot(rnd_param,std_RND*(ratio))
% plot(rnd_param,graphs.get_vertice_std('gaussian',rnd_param,signal_params.N));
o = fitoptions('poly3');
f = fit(rnd_param',std_RND,'poly3',o)
% plot(rnd_param,f.p1*rnd_param.^3+ f.p2*rnd_param.^2 + f.p3*rnd_param + f.p4);