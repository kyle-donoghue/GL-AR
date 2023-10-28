function coeffs = approxAR(y, P)
    mdl = varm(1,P);
    mdl.Constant = 0;
    for i = 1:size(y,2)
        estMdl = estimate(mdl,y(:,i));
        coeffs(:,i) = [1;-horzcat(estMdl.AR{:})';zeros(P-estMdl.P,1)];
    end
end