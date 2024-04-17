function [coeffs, y, x, A] = sineBA(n, P, delay, interval, N, Fs, freqs, noise)
    %n - node count
    %P - AR approx order
    %interval - samples per buffer
    %inov - stimuli sample count
    %N - total sample count (buffer count = N/interval)
    %Fs - sampling freqs as column vec
    %freqs - node stimuli frequencies
    %noise - noise variance
    x = createSine(n,N,Fs,freqs,noise);
    interval_length = floor(N/interval);
    A = adjacencyBA(n);
    if delay == 0
        D_AR = eye(n)+A;
    else
        mdl = createVARM(n, delay, A);
        D_AR = cat(3,mdl.AR{:});
        D_AR(:,:,2:end+1) = D_AR(:,:,1:end);
        D_AR(:,:,1) = eye(n);
    end
    y = zeros(interval_length,n,interval);
    coeffs = zeros(P+1,n,interval);
    for i = 1:interval
        % y(:,:,i) = filter(mdl,x(((i-1)*interval_length+1):(i*interval_length),:));
        y(:,:,i) = customFilter(D_AR,x(((i-1)*interval_length+1):(i*interval_length),:),n);
        coeffs(:,:,i) = approxAR(y(:,:,i),P);
    end
end