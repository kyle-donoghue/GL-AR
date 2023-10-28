clear;close all;
addpath("func")
addpath("toolbox")
addpath("toolbox/HalfVectorization")
addpath("sgwt_toolbox")
addpath("sgwt_toolbox/utils")
addpath("../../optimization")
addpath("../../creating_signals")
%% Generate a graph
nreplicate = 10; % repeat the same experiment (based on different graphs)
for ii = 1:nreplicate

param.N = 20;
% [A,XCoords, YCoords] = construct_graph(param.N,'gaussian',0.75,0.5);
[A,XCoords, YCoords] = construct_graph(param.N,'er',0.2);
% [A,XCoords, YCoords] = construct_graph(param.N,'pa',1);

%% Generate the graph Laplacian 
L_01 = full(sgwt_laplacian(A,'opt','raw'));
L_0 = L_01/trace(L_01)*param.N;

%% generate training signals
[V,D] = eig(full(L_0));
sigma = pinv(D);
mu = zeros(1,param.N);
num_of_signal = 100;
gftcoeff = mvnrnd(mu,sigma,num_of_signal);
X = V*gftcoeff';
X_noisy = X + 0.5*randn(size(X));

%% set parameters
param.max_iter = 50;
% alpha = 10.^[-1:-0.1:-3];
% beta = 10.^[0:-0.1:-2];
% lambda = 10.^[3:-0.05:0];
alpha = 10.^(-2);
beta = 10.^(-0.2);
lambda = 10.^1;

precision = zeros(length(alpha),length(beta));
recall = zeros(length(alpha),length(beta));
Fmeasure = zeros(length(alpha),length(beta));
NMI = zeros(length(alpha),length(beta));
num_of_edges = zeros(length(alpha),length(beta));

precision_logdet = zeros(1,length(lambda));
recall_logdet = zeros(1,length(lambda));
Fmeasure_logdet = zeros(1,length(lambda));
NMI_logdet = zeros(1,length(lambda));
num_of_edges_logdet = zeros(1,length(lambda));

%% main loop
for i = 1:length(alpha)
    for j = 1:length(beta)
        param.alpha = alpha(i);
        param.beta = beta(j);
        % GL-SigRep
        [L,Y,~] = graph_learning_gaussian(X_noisy,param);
        Lcell{i,j} = L;
        L(abs(L)<10^(-4))=0;
        [precision(i,j),recall(i,j),Fmeasure(i,j),NMI(i,j),num_of_edges(i,j)] = graph_learning_perf_eval(L_0,L);
    end
end
%adding new stuff below:
for i = 1:length(alpha)
    for j = 1:length(beta)
        param.alpha = alpha(i);
        param.beta = beta(j);
        % GL-SigRep
        B = approxAR(X_noisy',30);
        L_new = opt_L(B',param.alpha,param.beta);
        Lcell_new{i,j} = L_new;
        L_new(abs(L_new)<10^(-4))=0;
        [precision_new(i,j),recall_new(i,j),Fmeasure_new(i,j),NMI_new(i,j),num_of_edges_new(i,j)] = graph_learning_perf_eval(L_0,L_new);
    end
end
for k = 1:length(lambda)
    param.lambda = lambda(k);
    % GL-LogDet
    L_logdet = graph_learning_logdet_reglap(X_noisy,param);
    Lcell_logdet{k} = L_logdet;
    L_logdet(abs(L_logdet)<10^(-4))=0;
    [precision_logdet(k),recall_logdet(k),Fmeasure_logdet(k),NMI_logdet(k),num_of_edges_logdet(k)] = graph_learning_perf_eval(L_0,L_logdet);
end

%% performance
result1(:,:,ii) = precision;
result2(:,:,ii) = recall;
result3(:,:,ii) = Fmeasure;
result4(:,:,ii) = NMI;
result5(:,:,ii) = num_of_edges;

result1_new(:,:,ii) = precision_new;
result2_new(:,:,ii) = recall_new;
result3_new(:,:,ii) = Fmeasure_new;
result4_new(:,:,ii) = NMI_new;
result5_new(:,:,ii) = num_of_edges_new;

result1_logdet(ii,:) = precision_logdet;
result2_logdet(ii,:) = recall_logdet;
result3_logdet(ii,:) = Fmeasure_logdet;
result4_logdet(ii,:) = NMI_logdet;
result5_logdet(ii,:) = num_of_edges_logdet;

data{ii,1} = A;
data{ii,2} = L_0;
data{ii,3} = X;
data{ii,4} = X_noisy;

graph_original{ii} = L_0;
graph{ii} = Lcell;
graph_logdet{ii} = Lcell_logdet;

end