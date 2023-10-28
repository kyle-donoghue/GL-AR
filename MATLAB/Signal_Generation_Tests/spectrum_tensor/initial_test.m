clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
N = 4;
p = 2;
l = 64;
M = l*p;
x = zeros(N,M);
x(1,6:30) = sin(.25*(1:25));
% x(1,6:10) = 1;

[~,w] = freqz(ones(20,1),1,64,'whole');
H2 = 1+exp(-j*w)+exp(-j*2*w)+exp(-j*3*w);
H2 = (1-3*exp(-j*w))./(1-.4*exp(-j*w));
G = zeros(N,N,l);
X = zeros(N,p,l);

G(2,1,:) = H2;

for i = 1:p
    tmp = fft(x(:,(1+(i-1)*M/p):(i*M/p))',l);
    X(:,i,:) = tmp.';
end

% for i = 1:p
%     for j = 1:N
%         tmp = reshape(X(j,i,:).*G(j,i,:),[l 1]);
%         x_new(j,:,i) = ifft(tmp)';
%     end
% end
% x_new = [x_new(:,:,1) x_new(:,:,2)];

for i = 1:l
    X_new(:,:,i) = G(:,:,i)*X(:,:,i);
end

for i = 1:p
    tmp = permute(X_new(:,i,:),[3 1 2]);
    x_new(:,:,i) = real(ifft(tmp))';
end

x_new = real([x_new(:,:,1) x_new(:,:,2)]);

figure;
subplot(1,2,1)
signals.plot(x')
subplot(1,2,2)
signals.plot(x_new')

figure;hold;
for i = 1:100
    zero = randn(1);
    pole = .5*randn(1);
    H = (1-zero*exp(-j*w))./(1-pole*exp(-j*w));
    plot(10*log10(abs(H)));
end
figure
plot(abs(H))