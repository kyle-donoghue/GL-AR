function A = adjacencyER_random(n,p)
    if n == 1
        A = zeros(1,1);
        return;
    end
    A = rand(n,n) < p;
    A = double(triu(A,1));
    for i = 1:n
        for j = i:n
            if A(i,j) > 0
                A(i,j) = 2*rand(1)-1;
            end
        end
    end
    A = A + A';
    max_eig = max(abs(eig(A)));
    if max_eig ~= 0
        A = A/(1.1*max_eig);
    end
end