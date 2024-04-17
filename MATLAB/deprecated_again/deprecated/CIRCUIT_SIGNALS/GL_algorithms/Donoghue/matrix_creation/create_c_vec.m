function c = create_c_vec(N,Y)
    M = create_dup_matrix(N);
    Y2 = Y*Y';
    c = M'*Y2(:);
end