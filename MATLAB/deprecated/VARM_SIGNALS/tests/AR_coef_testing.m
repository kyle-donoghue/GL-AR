close all hidden;
load("x_data_for_ar.mat");

P = 10;

idp1 = ar(X(:,1)/norm(X(:,1)),P,'yw');
idp14 = ar(X(:,14)/norm(X(:,14)),P,'yw');
idp15 = ar(X(:,15)/norm(X(:,15)),P,'yw');

figure
plot(X(:,[1 14 15]))
figure;hold;
plot(idp1.A)
plot(idp14.A)
plot(idp15.A)

figure
freqz(1,idp1.A);
figure
freqz(1,idp14.A);
figure
freqz(1,idp15.A);

p_order = 2:2:45;
aic = zeros(length(p_order),1);
for i = 1:length(p_order)
    aic(i) = ar(X(:,15),p_order(i),'yw').Report.Fit.AIC;
end
figure
stem(p_order,aic)
figure;hold
f1 = abs(fft(X(:,1),128));
f14 = abs(fft(X(:,14),128));
f15 = abs(fft(X(:,15),128));
plot(f1/norm(f1))
plot(f14/norm(f14))
plot(f15/norm(f15))
