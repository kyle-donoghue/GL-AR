function A = adjacencyER(n,p)
    A = rand(n,n) < p;
    A = triu(A,1);
    A = A + A';