function x = createPulse(freq,Fs)
    N = floor(Fs/freq);
    t = 0:1/Fs:(N-1)/Fs;
    x = sin(2*pi*freq*t);
end