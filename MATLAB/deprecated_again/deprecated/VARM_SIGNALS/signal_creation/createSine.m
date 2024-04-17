function x = createSine(n,N,Fs,freqs,noise)
    t = 1/Fs:1/Fs:N/Fs;
    x = sin(2*pi*freqs*t)'+wgn(N,n,noise,'linear');
end