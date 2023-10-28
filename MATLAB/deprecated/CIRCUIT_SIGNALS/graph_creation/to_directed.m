function A = to_directed(A)
    n = size(A,1);
    for i = 1:n
        for j = i+1:n
            if A(i,j) > 0
                r = rand(1) > .5;
                if r
                    A(i,j) = 0;
                else
                    A(j,i) = 0;
                end
            end
        end
    end
end