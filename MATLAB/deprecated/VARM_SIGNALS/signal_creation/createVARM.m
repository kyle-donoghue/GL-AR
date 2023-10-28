function mdl = createVARM(n, P, A)
    mdl = varm(n,P);
    mdl.Constant = zeros(n,1);
    for i = 1:P-1
        mdl.AR(i) = {zeros(n,n)};
    end
    mdl.AR(P) = {A}; %set P-delay to adjacency
    mdl.Covariance = eye(n);
end