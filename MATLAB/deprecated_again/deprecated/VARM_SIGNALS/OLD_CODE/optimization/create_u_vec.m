function u = create_u_vec(N,i,j)
    u = zeros(N*(N+1)/2,1);
    u((j-1)*N+i-j*(j-1)/2) = 1;
end