clc;clear;close all hidden;
N = 20;
[A,XCoords, YCoords] = construct_graph(N,'er',0.2);
L = adjacency_to_laplacian(A);
t1 = trace(L)
s1 = sum(eig(A).^2)

L_1 = L/trace(L)*(N-1);
A_1 = laplacian_to_adjacency(L_1);

t = trace(L_1)
s = sum((eig(A_1)).^2)
s2 = trace(L)/(N-1)*sum(eig(A_1).^2)

% s2 should equal N-1 to max that
% which means trace(L/alpha) = alpha*sum(eig(A).^2)
% we want sum(eig(A).^2) < N or = N-1
% so trace(L/alpha)/alpha = N-1
% or trace(L/alpha) = (N-1)*alpha