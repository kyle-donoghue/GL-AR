function A = adjacencyRND(n,r)
    A = zeros(n,n);
    if n == 1
        return
    end
    cords = rand(2,n);
    for i = 1:n
        for j = 1:i-1
            if norm(cords(:,i)-cords(:,j)) < r
                A(i,j) = 1;
                A(j,i) = 1;
            end
        end
    end
    max_eig = max(abs(eig(A)));
    if max_eig ~= 0
        A = A/(1.1*max_eig);
    end
end