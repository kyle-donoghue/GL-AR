function [L_dist, A_dist] = create_temperature_groundtruth(DIST,sigma,alpha)
    A_dist = zeros(length(DIST));
    for m = 1:length(DIST)
        for n = 1:length(DIST)
            if m==n
                continue;
            end
            A_dist(m,n) = exp(-(DIST(m,n)/sigma)^2);
            denom1 = 0;
            denom2 = 0;
            [sort_row,ind1] = sort(DIST(m,:),'ascend');
            [sort_col,ind2] = sort(DIST(:,n),'ascend');
            for i = 1:alpha
                denom1 = denom1+exp(-(DIST(m,ind1(i))/sigma)^2);
                denom2 = denom2+exp(-(DIST(ind2(i),n)/sigma)^2);
            end
            A_dist(m,n) = A_dist(m,n)/sqrt((denom1*denom2));
        end
    end
    L_dist = graphs.to_laplacian(A_dist);
    L_dist = GL.threshold(L_dist,1e-2);
    A_dist = graphs.to_adjacency(L_dist);
end