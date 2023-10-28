clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_approximation/"));
load("../signal_creation/n20.mat")
figure
matrixPlot(A)
figure
signalPlot(out')