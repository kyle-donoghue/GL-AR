%%
clc;clear;close all hidden;
%%
% 1=sqrt(1^2*gauss+gamma^2*exp) ??
% pdf = sqrt(1/sqrt(1+gamma^2)*gauss + gamma/sqrt(1+gamma^2)*exp) ?
% nope! its this:
% f(x) = 1/(1+gamma)*f_gauss(x)+gamma/(1+gamma)*f_exp(x)
% i think the above is correct, but i think lambda for exp and variance for gauss should scale with gamma too
% for ex, lambda goes down as gamma increases, and variance goes down as gamma decreases
% sigma^2 = sigma_o^2*(1+gamma)
% lamba = lambda_o/(1+gamma)

lambda_0 = 10;
sigmasqr_0 = .0001;
mu = .05;
x =0:.001:1;

figure;
for gamma = logspace(-2,2,20)
    f_E = lambda_0*(1+gamma)*exp(-lambda_0*(1+gamma).*x);
    f_G = 1/sqrt(2*pi*sigmasqr_0*(1+gamma))*exp(-(x-mu).^2/(2*sigmasqr_0*(1+gamma)));
    f = 1/(1+gamma)^2*(f_G+gamma^2*f_E);
    plot(x,f)
    pause(.5)
end