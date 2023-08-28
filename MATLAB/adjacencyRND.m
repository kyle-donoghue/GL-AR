function A = adjacencyRND(n,r)
    A = zeros(n,n);
    cords = rand(2,n);
    for i = 1:n
        for j = 1:i-1
            if norm(cords(:,i)-cords(:,j)) < r
                A(i,j) = 1;
                A(j,i) = 1;
            end
        end
    end
end