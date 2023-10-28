function [W,sigmasquare] = optimize_laplacian_logdet_reglap(Y,lambda)

%% Laplacian constraints
[N,m] = size(Y);
%% optimization
cvx_begin

% cvx_solver mosek

variable Delta1(N,N) semidefinite
variable W(N,N) symmetric
variable sigmasquare

minimize trace(1/m*Delta1*Y*Y') - log_det(Delta1) + 1/m*lambda*norm(reshape(W,N^2,1),1)

subject to
    Delta1 == diag(sum(W,2)) - W + eye(N)*sigmasquare % sigmasquare is the inverse of sigma^2 in the paper
    diag(W) == 0
    reshape(W,N^2,1) >= 0
    sigmasquare >= 10^(-4) % because strict inequality is discouraged

cvx_end