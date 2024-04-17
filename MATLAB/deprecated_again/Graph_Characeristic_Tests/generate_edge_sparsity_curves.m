% clc;clear;close all hidden;
addpath(genpath('../GL_classes/'));
load("edge_sparsity.mat");
%% ER Curve Fit
figure;
subplot(1,2,1); hold on;
plot(er_param,sparsity_ER);
plot(er_param,graphs.get_sparsity('er',er_param));
subplot(1,2,2); hold on;
plot(er_param,edge_spread_ER)
plot(er_param,graphs.get_edge_spread('er',er_param,20));

%% BA Curve Fit
figure;
subplot(1,2,1); hold on;
plot(ba_param,sparsity_BA);
% o = fitoptions('poly2');
% o.Lower = [-Inf -Inf 0];
% o.Upper = [Inf Inf 0];
% f = fit(ba_param',sparsity_BA,'poly2',o)
% plot(ba_param,f.p1*ba_param.^2 + f.p2*ba_param + f.p3);
plot(ba_param,graphs.get_sparsity('pa',ba_param));
subplot(1,2,2); hold on;
plot(ba_param,edge_spread_BA)
% f = fit(ba_param',edge_spread_BA/20,'poly2')
% plot(ba_param,20*(f.p1*ba_param.^2 + f.p2*ba_param + f.p3));
plot(ba_param,graphs.get_edge_spread('pa',ba_param,10));

%% RND Curve Fit
figure;
subplot(1,2,1); hold on;
plot(rnd_param,sparsity_RND);
% o = fitoptions('poly3');
% o.Lower = [-Inf -Inf -Inf 0];
% o.Upper = [Inf Inf Inf 0];
% f = fit(rnd_param',sparsity_RND,'poly3',o)
% plot(rnd_param,f.p1*rnd_param.^3+ f.p2*rnd_param.^2 + f.p3*rnd_param + f.p4);
plot(rnd_param,graphs.get_sparsity('gaussian',rnd_param));
subplot(1,2,2); hold on;
plot(rnd_param,edge_spread_RND)
% f = fit(rnd_param',edge_spread_RND/20,'poly3',o)
% plot(rnd_param,(f.p1*rnd_param.^3+ f.p2*rnd_param.^2 + f.p3*rnd_param + f.p4)*20);
plot(rnd_param,graphs.get_edge_spread('gaussian',rnd_param,20));
