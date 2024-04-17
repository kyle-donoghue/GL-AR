function c = create_c_vec(N,Y,alpha)
    M = create_dup_matrix(N);
    Y2 = Y*Y';
    c = alpha*M'*Y2(:);
end