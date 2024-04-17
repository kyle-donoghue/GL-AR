function A = adjacencyER(n,p)
    if n == 1
        A = zeros(1,1);
        return;
    end
    A = rand(n,n) < p;
    A = triu(A,1);
    A = A + A';
    max_eig = max(abs(eig(A)));
    if max_eig ~= 0
        A = A/(1.1*max_eig);
    end
end