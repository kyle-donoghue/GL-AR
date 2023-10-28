function A = laplacian_to_adjacency(L)
    A = diag(diag(L))-L;
end