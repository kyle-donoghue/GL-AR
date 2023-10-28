clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
rng('default')
%%
signal_params.N = 8;
signal_params.raw = 1;

signal_params.SNR = 2e0;
signal_params.Fs = 100;
signal_params.min_freq = 2;
signal_params.max_freq = 45;

signal_params.intervals = signal_params.N;
signal_params.interval_length = 128;
signal_params.M = signal_params.intervals*signal_params.interval_length;

signal_params.order = 1;
signal_params.zero_variance = 1;
signal_params.zero_mean = 1;
signal_params.pole_variance = .25;
signal_params.pole_mean = .5;

%%
AR_params.N = signal_params.N;
AR_params.P = 10;
AR_params.gamma = 100;
AR_params.threshold = 0;
%% 
[L_0,A] = graphs.create(signal_params,'er',.2);
A_0 = graphs.to_directed(A);
G = graphs.createGraphTensor(signal_params,A_0);
x = randn(signal_params.N,signal_params.M);
x = signals.create_raw_sine(signal_params);
X = signals.createTensorSpectrum(signal_params,x);
Y = signals.filterTensorSpectrum(signal_params,X,G);
y = signals.inverseTensorSpectrum(signal_params,Y);
y_noisy = signals.add_noise(signal_params,y);

figure;
subplot(2,2,1);plot(x(2,1:signal_params.interval_length)');
subplot(2,2,2);plot(y(2,1:signal_params.interval_length)');
subplot(2,2,3);plot(x(3,1:signal_params.interval_length)');
subplot(2,2,4);plot(y(3,1:signal_params.interval_length)');
figure;
subplot(1,2,1);plot(abs(fft(x(2,1:signal_params.interval_length)')));
subplot(1,2,2);plot(abs(fft(y(2,1:signal_params.interval_length)')));

figure;
subplot(1,3,1);signals.plot(signals.approxAR(x(:,1:signal_params.interval_length),AR_params.P)')
subplot(1,3,2);signals.plot(signals.approxAR(y_noisy(:,1:signal_params.interval_length),AR_params.P)')
B = zeros(signal_params.N,AR_params.P,signal_params.intervals);
for i = 1:signal_params.intervals
    B(:,:,i) = signals.approxAR(y_noisy(:,(1+(i-1)*signal_params.interval_length):(i*signal_params.interval_length)),AR_params.P);
end
B = mean(B,3);
subplot(1,3,3);signals.plot(B')

[fits,aic,p] = signals.test_fit(y_noisy(1,1:signal_params.interval_length));
[fits2,aic2,p2] = signals.test_fit(x(1,1:signal_params.interval_length));

figure;
signals.plot([x(1,1:signal_params.interval_length)' y(1,1:signal_params.interval_length)']);

figure;
signals.plot([fits fits2]);

%%
L = zeros(signal_params.N,signal_params.N,signal_params.intervals);
for i = 1:signal_params.intervals
    L(:,:,i) = GL.AR(y_noisy(:,(1+(i-1)*signal_params.interval_length):(i*signal_params.interval_length)),AR_params);
end
L = mean(L,3);
%%
p = AR_params;
Q = create_Q_matrix(p.N,1);

c = create_c_vec(p.N,B)*p.gamma;
A = create_constraint_matrix(p.N);
b = create_constraint_vec(p.N,p.N);

% options = optimset('Display', 'off');
% [phi,f] = quadprog(Q,c,A(p.N+2:end,:),b(p.N+2:end),A(1:p.N+1,:),b(1:p.N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...

lowerbound = -Inf(size(b));
lowerbound(1:p.N+1) = b(1:p.N+1);
upperbound = b;
m = osqp;
m.setup(Q,c,A,lowerbound,upperbound,'verbose',false);
results = m.solve();
phi = results.x;

l = create_dup_matrix(p.N)*phi;
L2 = convert_to_matrix(l);
%%
figure;
graphs.vcompare(L_0,L2,'L_0','L2');