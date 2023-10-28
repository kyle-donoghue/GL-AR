clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
rng("default")
%% Signal Configuration
signal_params.N = 20;
signal_params.M = 1e3;
signal_params.active = 20;
signal_params.SNR = 2;
signal_params.additional_diag = 2.5; %shunt resistors
signal_params.trace_normalization = signal_params.N;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;
signal_params.isComplex = 1;

%% GL_AR Configuration
AR_params.N = signal_params.N;
AR_params.P = 20;
AR_params.gamma = 20;%2.6;
AR_params.threshold = 1e-3;%.085;

%% Dong Configuration
dong_params.N = signal_params.N;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.beta = 10.^(-3);
dong_params.threshold = .07;

%% Graph Creation
[L_0,A_0] = graphs.create(signal_params,'er',0.2);
figure;graphs.plot(A_0,getVarName(A_0));

%% Signal Creation
[V,I, active_indices,freqs] = signals.circuit_sine(signal_params,L_0);
% V = signals.z_score(V);

%% GL_AR
L_AR = GL.AR(V,AR_params);
figure;graphs.plot(L_AR,getVarName(L_AR))
L_AR = GL.threshold(L_AR,AR_params.threshold);

%% GL_Dong
[L_dong,Y,~] = GL.dong(V,dong_params);
figure;graphs.plot(L_dong,getVarName(L_dong))
L_dong = GL.threshold(L_dong,dong_params.threshold);

%% Comparison
[p_AR,r_AR,f_AR,nmi_AR,edges_AR] = graphs.performance(L_0,L_AR);
[p_dong,r_dong,f_dong,nmi_dong,edges_dong] = graphs.performance(L_0,L_dong);

A_AR = graphs.to_adjacency(L_AR);
A_dong = graphs.to_adjacency(L_dong);

figure;
graphs.vcompare(abs(A_0),abs(A_AR),getVarName(A_0),getVarName(A_AR));
figure;
graphs.vcompare(abs(A_0),abs(A_dong),getVarName(A_0),getVarName(A_dong));
figure;
graphs.vcompare(abs(A_AR),abs(A_dong),getVarName(A_AR),getVarName(A_dong));