clc;clear;close all hidden;
load('large_maximum_data_results_n20_n10.mat');
load('large_maximum_data_results_n40.mat');
%%
N=20;
sparse_F = F_RND_20;
sparsities = S_RND_20;
thresholds = 0:.0025:2;
edge_spreads = graphs.get_edge_spread('gaussian',rnd_param,N);
xaxis = edge_spreads/N;
max_thresholds = zeros(length(sparsities),1);
z = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [z(i), ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    max_thresholds(i) = thresholds(col);
end
max_thresholds = max_thresholds*N;
%%
max_thresholds(8:end) = .25;
%%
figure;hold on;
plot(max_thresholds);
%%
f = fit(xaxis',max_thresholds-.25,'exp2')
y_hat = f.a*exp(f.b*xaxis) + f.c*exp(f.d*xaxis)+.25;
plot(y_hat)
%%
N=40;
sparse_F = F_RND_40;
sparsities = S_RND_40;
thresholds = 0:.0025:2;
edge_spreads = graphs.get_edge_spread('gaussian',rnd_param,N);
xaxis = edge_spreads/N;
max_thresholds = zeros(length(sparsities),1);
z = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [z(i), ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    max_thresholds(i) = thresholds(col);
end
max_thresholds = max_thresholds*N;
max_thresholds(7:end) = .25;

%%
plot(max_thresholds);
o = fitoptions('exp2');
o.Lower = [f.a f.b -Inf -Inf];
o.Upper = [f.a f.b Inf Inf];

f2 = fit(xaxis',max_thresholds-.25,'exp2',o)
y_hat = f2.a*exp(f2.b*xaxis) + f2.c*exp(f2.d*xaxis)+.25;
plot(y_hat)
y_hat = f2.a*exp(f2.b*xaxis) + f.c/2*exp(f.d*sqrt(2)*xaxis)+.25;
plot(y_hat,'--')
