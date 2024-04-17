function [coeffs, y, x, A] = sineER_random(n, P, p, delay, interval, N, Fs, freqs, noise)
    %n - node count
    %P - AR approx order
    %interval - samples per buffer
    %inov - stimuli sample count
    %N - total sample count (buffer count = N/interval)
    %Fs - sampling freqs as column vec
    %freqs - node stimuli frequencies
    %noise - noise variance
    A = adjacencyER_random(n,p);
    mdl = createVARM(n, delay, A);
    x = createSine(n,N,Fs,freqs,noise);
    interval_length = floor(N/interval);
    D_AR = cat(3,mdl.AR{:});
    D_AR(:,:,2:end+1) = D_AR(:,:,1:end);
    D_AR(:,:,1) = eye(n);
    for i = 1:interval
        % y(:,:,i) = filter(mdl,x(((i-1)*interval_length+1):(i*interval_length),:));
        y(:,:,i) = customFilter(D_AR,x(((i-1)*interval_length+1):(i*interval_length),:),n);
        coeffs(:,:,i) = approxAR(y(:,:,i),P);
    end
end