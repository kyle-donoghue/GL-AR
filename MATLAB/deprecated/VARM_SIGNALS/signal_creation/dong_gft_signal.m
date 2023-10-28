function X = dong_gft_signal(L,p)
    l = p.l;
    N = p.N;
    noise = p.noise;
    [V,D] = eig(full(L));
    sigma = pinv(D);
    mu = zeros(1,N);
    gftcoeff = mvnrnd(mu,sigma,l);
    x = V*gftcoeff';
    X = x + noise*randn(size(x));
end