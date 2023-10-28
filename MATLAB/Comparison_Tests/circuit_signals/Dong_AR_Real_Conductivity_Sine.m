clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%% Signal Configuration
signal_params.N = 20;
signal_params.M = 1e3;
signal_params.active = 20;
signal_params.SNR = 2;
signal_params.additional_diag = .15; %shunt resistors
signal_params.trace_normalization = signal_params.N;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;
signal_params.isComplex = 0;

%% GL_AR Configuration
AR_params.N = signal_params.N;
AR_params.P = 20;
AR_params.gamma = 10;
AR_params.threshold = .1;

%% Dong Configuration
dong_params.N = signal_params.N;
dong_params.max_iter = 50;
dong_params.alpha = 10.^(-2);
dong_params.beta = 10.^(.7);
dong_params.threshold = 10^-1;

%% Graph Creation
L_0 = graphs.create(signal_params,'er',0.2);

%% Signal Creation
V = signals.circuit_sine(signal_params,L_0);

%% GL_AR
L_AR = GL.AR(V,AR_params);
L_AR = GL.threshold(L_AR,AR_params.threshold);

%% GL_Dong
[L_dong,Y,~] = GL.dong(V,dong_params);
L_dong = GL.threshold(L_dong,dong_params.threshold);

%% Comparison
[p_AR,r_AR,f_AR,nmi_AR,edges_AR] = graphs.performance(L_0,L_AR);
[p_dong,r_dong,f_dong,nmi_dong,edges_dong] = graphs.performance(L_0,L_dong);

A_AR = graphs.to_adjacency(L_AR);
A_dong = graphs.to_adjacency(L_dong);
A_0 = graphs.to_adjacency(L_0);

figure;
graphs.vcompare(A_0,A_AR,'Groundtruth','AR-Recovered');
figure;
graphs.vcompare(A_0,A_dong,'Groundtruth','Dong-Recovered');
figure;
graphs.vcompare(A_AR,A_dong,'AR-Recovered','Dong-Recovered');