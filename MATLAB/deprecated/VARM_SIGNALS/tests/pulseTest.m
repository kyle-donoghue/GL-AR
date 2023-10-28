freqs = generate_frequencies2(global_params.N,2,25);
noise = 0;
N = 20;
[A,XCoords, YCoords] = construct_graph(N,'er',0.2);
A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A))));
A_0 = triu(A_0);
l = 1e3;
Fs = 100;
delay = 10;
minSep = 100;
maxSep = 500;
x = createFullSinePulse(N,l,Fs,freqs,noise,minSep,maxSep);
D_AR = zeros(N,N,delay+1);
D_AR(:,:,1) = eye(N);
D_AR(:,:,delay+1) = A_0;
y = customFilter(D_AR,x,N);
plot(y)