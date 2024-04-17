function y = varm_pulse_signal(A,params,freqs)
    %n - node count
    %P - AR approx order
    %interval - samples per buffer
    %inov - stimuli sample count
    %N - total sample count (buffer count = N/interval)
    %Fs - sampling freqs as column vec
    %freqs - node stimuli frequencies
    %noise - noise variance
    noise = params.noise;
    %freqs = params.freqs;
    N = params.N;
    l = params.l;
    Fs = params.Fs;
    delay = params.delay;
    minSep = params.minSep;
    maxSep = params.maxSep;
    x = createFullSinePulse(N,l,Fs,freqs,noise,minSep,maxSep);
    D_AR = zeros(N,N,delay+1);
    D_AR(:,:,1) = eye(N);
    D_AR(:,:,delay+1) = A;
    y = customFilter(D_AR,x,N);
    y = y';
end