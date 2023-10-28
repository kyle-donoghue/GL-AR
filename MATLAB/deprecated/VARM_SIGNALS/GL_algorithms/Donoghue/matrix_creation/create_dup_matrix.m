function M = create_dup_matrix(N)
    M = zeros(N*(N+1)/2,N^2);
    for j = 1:N
        for i = j:N
            T = create_T_matrix(N,i,j);
            M = M+create_u_vec(N,i,j)*T(:)';
        end
    end
    M = M';
end