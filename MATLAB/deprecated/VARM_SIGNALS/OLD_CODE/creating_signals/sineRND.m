function [coeffs, y, x, A] = sineBA(n, P, r, delay, interval, N, Fs, freqs, noise)
    %n - node count
    %P - AR approx order
    %interval - samples per buffer
    %inov - stimuli sample count
    %N - total sample count (buffer count = N/interval)
    %Fs - sampling freqs as column vec
    %freqs - node stimuli frequencies
    %noise - noise variance
    A = adjacencyRND(n,r);
    %A = [0, 1, 0, 0; 0, 0, 1, 0;0, 0, 0, 0;0, 0, 0, 0];
    mdl = createVARM(n, delay, A);
    x = createSine(n,N,Fs,freqs,noise);
    y = filter(mdl,x);
    coeffs = approxAR(y,P);
end