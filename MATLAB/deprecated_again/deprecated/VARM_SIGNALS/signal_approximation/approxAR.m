function coeffs = approxAR(y, P)
    % mdl = varm(1,P);
    % mdl.Constant = 0;
    coeffs = zeros(P+1,size(y,2));
    for i = 1:size(y,2)
        % estMdl = estimate(mdl,y(:,i));
        % coeffs(:,i) = [1;-horzcat(estMdl.AR{:})';zeros(P-estMdl.P,1)]; 
        coeffs(:,i) = ar(y(:,i),P,'yw').A; % switched to this for faster speed
        % coeffs(:,i) = ar(y(:,i),P,'ls').A; % switched to this for faster speed
    end
    coeffs = coeffs(2:end,:);
end