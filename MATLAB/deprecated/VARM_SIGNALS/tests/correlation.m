clc;clear;close all hidden;

% for i = 1:100
    x = randn(1,1e2);
    y = randn(1,1e2);
    x_a = ar(x,30,'ls').A;
    y_a = ar(y,30,'ls').A;
    x_a = x_a(2:end);
    y_a = y_a(2:end);
    e = x*y';
    e_a = x_a*y_a';
% end

% e = e./100;
% e_a = e_a./100;