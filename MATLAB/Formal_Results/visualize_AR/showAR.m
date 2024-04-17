clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
x = zeros(2,100);
x(1,20:40) = 1;
x(2,65:85) = 1;
B = signals.approxAR(x,40);
plot(B(1,:),'LineWidth',5)
hold on
plot(B(2,:)-1.5,'LineWidth',5)
axis off


hold off
plot(x(1,:),'LineWidth',5)
hold on
plot(x(2,:)-1.5,'LineWidth',5)
axis off