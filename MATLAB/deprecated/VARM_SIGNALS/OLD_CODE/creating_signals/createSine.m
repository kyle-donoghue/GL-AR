function x = createSine(n,N,Fs,freqs,noise)
    t = 1/Fs:1/Fs:N/Fs;
    x = sin(2*pi*freqs*t)'+wgn(N,n,noise,'linear');
    % x(1:innov,:) = sin(2*pi*freqs*t)'+wgn(innov,n,noise,'linear');
    %x(innov+1:N,:) = wgn(N-innov,n,noise,'linear');
end