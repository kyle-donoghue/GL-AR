function Q = create_Q_matrix(N,gamma)
    M = create_dup_matrix(N);
    Q = 2*gamma*M'*M;
end