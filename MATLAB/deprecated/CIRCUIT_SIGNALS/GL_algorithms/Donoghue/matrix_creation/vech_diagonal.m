function d = vech_diagonal(N)
    d = zeros(N,1);
    for j = 0:N-1
        d(j+1) = 1+j*N-(j-1)*j/2;
    end
end