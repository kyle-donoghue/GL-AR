function [L,iterations] = full_optimization(X,alpha,beta,tolerance)
    Y = X;
    iterations = 1;
    tol = inf;
    L = opt_L(Y,alpha,beta);
    Y = opt_Y(X,L,alpha);
    while tol > tolerance
        L_i = opt_L(Y,alpha,beta);
        Y = opt_Y(X,L_i,alpha);
        tol = max(abs(L(:)-L_i(:)));
        L = L_i;
        iterations = iterations+1;
    end
end