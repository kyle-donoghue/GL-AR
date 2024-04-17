function L = opt_L(Y, alpha, beta)
    N = size(Y,1);
    Q = create_Q_matrix(N,beta);
    c = create_c_vec(N,Y,alpha);
    A = create_constraint_matrix(N);
    b = create_constraint_vec(N);
    model.Q = sparse(Q);
    model.A = sparse(A);
    model.obj = c;
    model.rhs = b;
    model.sense = [char('='*ones(N+1,1));char('<'*ones(N*(N-1)/2,1))];
    model.lb = -inf*ones(N*(N+1)/2,1);
    params.outputflag = 0;
    results = gurobi(model,params);
    phi = results.x;
    l = create_dup_matrix(N)*phi;
    L = convert_to_matrix(l);
end