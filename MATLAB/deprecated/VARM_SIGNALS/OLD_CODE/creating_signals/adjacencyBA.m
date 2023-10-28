function A = adjacencyBA(n)
    if n == 1
        A = zeros(1,1);
        return
    end
    m = 1;
    %if (Nodes < pos) || (mlinks > pos) || (pos < 1) || (size(size(seed)) ~= 2) || (mlinks < 1) || (seed ~= seed') || (sum(diag(seed)) ~= 0)
    %    error('invalid parameter value(s)');
    %end
    %if mlinks > 5 || Nodes > 15000 || pos > 15000
    %    warning('Abnormally large value(s) may cause long processing time');
    %end
    A = zeros(n, n);
    A(1:2,1:2) = [0 1;1 0];
    sumlinks = sum(sum(A));
    pos = 2;
    while pos < n
        pos = pos + 1;
        linkage = 0;
        while linkage ~= m
            rnode = ceil(rand * pos);
            deg = sum(A(:,rnode)) * 2;
            rlink = rand * 1;
            if rlink < deg / sumlinks && A(pos,rnode) ~= 1 && A(rnode,pos) ~= 1
                A(pos,rnode) = 1;
                A(rnode,pos) = 1;
                linkage = linkage + 1;
                sumlinks = sumlinks + 2;
            end
        end
    end
    max_eig = max(abs(eig(A)));
    if max_eig ~= 0
        A = A/(1.1*max_eig);
    end
end
