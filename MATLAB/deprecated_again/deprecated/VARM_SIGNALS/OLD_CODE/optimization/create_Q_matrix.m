function Q = create_Q_matrix(N,beta)
    M = create_dup_matrix(N);
    Q = 2*beta*M'*M;
end