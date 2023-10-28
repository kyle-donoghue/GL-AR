function A =  create_constraint_matrix2(N)
    A = zeros((2*N+1)+N*(N-1)/2,N*(N+1)/2);
    % tr(L)=N
    diag_indices = zeros(N,1);
    for j = 0:N-1 % add diagonal elements based off vech indexing
        diag_indices(j+1) = 1+j*N-(j-1)*j/2;
    end
    for i = 1:N % make diagonal elements proportional to signal energy
        A(N+1+i,diag_indices(i)) = 1;
    end
    A(N+1,diag_indices) = -1; 
    % L*1=0
    for j = 1:N
        A(j,diag_indices(j):diag_indices(j)+(N-j)) = 1; % add elements below diagonal based off vech indexing
        left_indices = zeros(j-1,1);
        for jj = 1:j-1 % add elements left of diagonal based off vech indexing
            left_indices(jj) = diag_indices(j)-(N-j+1);
            for k = j-jj+2:j
                left_indices(jj) = left_indices(jj) - (N+2-k);
            end
        end
        A(j,left_indices) = 1;
    end
    % L_ij < 0
    i = 2*N+2;
    for k = 1:N
        for j = diag_indices(k)+1:diag_indices(k)+(N-k) % add off diagonal elements individually based off vech indexing
            A(i,j) = 1;
            i = i+1;
        end
    end
end