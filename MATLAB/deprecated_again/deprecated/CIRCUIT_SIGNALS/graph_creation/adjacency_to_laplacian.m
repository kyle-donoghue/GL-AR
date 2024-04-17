function L = adjacency_to_laplacian(A)
    L = diag(sum(A))-A;
end