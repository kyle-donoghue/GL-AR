function y = customFilter(D,x,m)
    p = length(D(1,1,:))-1;
    n = length(x);
    y = zeros(n,m);
    for i = 1:p
        y(i,:) = x(i,:);
        for k = 2:i
            y(i,:) = y(i,:) + (D(:,:,k)*y(i-k+1,:)')';
        end
    end
    for i = p+1:n
        y(i,:) = x(i,:);
        for k = 2:p+1
            y(i,:) = y(i,:) + (D(:,:,k)*y(i-k+1,:)')';
        end
        % y(i) = -D(2:end)'*flipud(y(i-p:i-1))+x(i);
    end
end