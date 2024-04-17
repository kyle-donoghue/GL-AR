function X = convert_to_matrix(x)
    X = zeros(sqrt(length(x)));
    for i = 1:length(X)
        for j = 1:length(X)
            X(j,i) = x((i-1)*length(X)+j);
        end
    end
end