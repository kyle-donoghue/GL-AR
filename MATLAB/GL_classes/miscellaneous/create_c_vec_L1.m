function c = create_c_vec_L1(N,Y,gamma)
    M = create_dup_matrix(N);
    Y2 = Y*Y';
    y = Y2(:);
    spars = ones(size(y));
    for i = 1:length(spars)
        if mod(i-1,N) == 0
            spars(i) = 0;
        end
    end
    c = M'*(gamma*y(:)+spars);
end