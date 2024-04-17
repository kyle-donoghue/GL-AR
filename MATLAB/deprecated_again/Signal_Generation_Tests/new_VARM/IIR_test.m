clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
p.pole_variance = .1;
p.pole_mean = .7;
p.order = 1;
%%
figure;hold;
for i = 1:20
    a = signals.create_IIR(p);
    plot(abs(freqz(1,a)));
end

figure;
zplane(1,a)