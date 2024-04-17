function T = create_T_matrix(N,i,j)
    T = sparse(N,N);
    T(i,j) = 1;
    T(j,i) = 1;
end