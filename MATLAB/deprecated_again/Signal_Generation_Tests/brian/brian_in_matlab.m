clc;clear;clear classes;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
% out = pyrunfile("create_brian.py","out",ML_A=A);
% V = double(out)';
% plot(V)
tic

A = zeros(4);
A(1,3) = 1;
dur = 80;
stimulus = [0, 2, .5, .04; 1, 2, .5, .04; 0, 10, .5, .15; 1, 25, .5, .15; 0, 35, .5, .15; 1, 45, .5, .15];
% parfor k = 1:2
    python.imports;
    [t,v,i] = python.simulate_HH_network(A,stimulus,dur);
% end
plot(t,v')
toc