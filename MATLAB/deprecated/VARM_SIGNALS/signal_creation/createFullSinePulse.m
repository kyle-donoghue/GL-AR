function x = createFullSinePulse(n,N,Fs,freqs,noise,minSeparation,maxSeparation)
    x = zeros(N,n);
    for i = 1:n
        inSignal = 1;
        f = freqs(i);
        indices = [randi(maxSeparation)];
        while inSignal
           indices(end+1) = indices(end)+floor(Fs/f)+minSeparation+randi(maxSeparation);
           inSignal = (indices(end)+floor(Fs/f)) < N;
        end
        indices = indices(1:end-1);
        for k = 1:length(indices)
            x_tmp = createPulse(f,Fs);
            x(indices(k):indices(k)+length(x_tmp)-1,i) = x_tmp;
        end
    end
    x = x+wgn(N,n,noise,'linear');
end