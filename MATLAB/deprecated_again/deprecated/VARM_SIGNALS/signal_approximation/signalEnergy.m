function [e,e_i] = signalEnergy(X)
    e = 0;
    e_i = zeros(size(X,1),1);
    for i = 1:size(X,1)
        e_i(i) = norm(X(i,:))^2;
        e = e+e_i(i);
    end
end