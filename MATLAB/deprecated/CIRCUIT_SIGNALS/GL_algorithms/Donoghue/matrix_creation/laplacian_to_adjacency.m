function A = laplacian_to_adjacency(L)
    A = zeros(size(L));
    for i = 1:size(L,3)
        A(:,:,i) = diag(diag(L(:,:,i)))-L(:,:,i);
    end
end