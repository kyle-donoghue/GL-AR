function b = create_constraint_vec(N,trace_constraint)
    b = sparse((N+1)+N*(N-1)/2,1);
    b(N+1) = trace_constraint; %switched from N temporarily
end