function Y = opt_Y(X,L,alpha)
    Y = (eye(size(L))+alpha*L)^-1*X;
end