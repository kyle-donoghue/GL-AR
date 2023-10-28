function b = create_constraint_vec(N)
    b = zeros((N+1)+N*(N-1)/2,1);
    b(N+1) = N;
end