function b = create_constraint_vec2(N,sig_energy,e)
    b = zeros((2*N+1)+N*(N-1)/2,1);
    b(N+1) = -min(sig_energy)*N;
    b(N+2:2*N+1) = sig_energy; %switched from N temporarily
end