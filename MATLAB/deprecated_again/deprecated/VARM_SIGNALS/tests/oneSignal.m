clc;clear;close all hidden;
addpath(genpath("../graph_creation/"));
addpath(genpath("../graph_accuracy/"));
addpath(genpath("../signal_creation/"));
addpath(genpath("../signal_approximation/"));
addpath(genpath("../GL_algorithms/"));
addpath(genpath("../progressBar/"));
addpath(genpath("../Comparison_Tests/"));
rng("default")
global_params.N = 20; % node count
global_params.l = 1e3; % signal length
global_params.noise = .000005; %noise power
global_params.trials = 50;
global_params.Fs = 100;
global_params.delay = 10;
global_params.minSep = 100;
global_params.maxSep = 500;

donoghue_params = global_params;
donoghue_params.P = 20; % AR approx order
donoghue_params.gamma = .002;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.threshold_avg = 0.2;
donoghue_params.threshold = 1;

[A,XCoords, YCoords] = construct_graph(global_params.N,'er',.2);
A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
L_0 = adjacency_to_laplacian(A_0);
A_1 = to_directed(A_0);
freqs = generate_frequencies2(global_params.N,2,25);
% freqs = 10.2:.2:14;

x = zeros(global_params.l,global_params.N);
for i = 1:global_params.N
    inSignal = 1;
    f = freqs(i);
    indices = [randi(global_params.maxSep)];
    while inSignal
       indices(end+1) = indices(end)+floor(global_params.Fs/f)+global_params.minSep+randi(global_params.maxSep);
       inSignal = (indices(end)+floor(global_params.Fs/f)) < global_params.l;
    end
    indices = indices(1:end-1);
    for k = 1:length(indices)
        x_tmp = createPulse(f,global_params.Fs);
        x(indices(k):indices(k)+length(x_tmp)-1,i) = x_tmp;
    end
end
x = x+wgn(global_params.l,global_params.N,global_params.noise,'linear');

D_AR = zeros(global_params.N,global_params.N,global_params.delay+1);
D_AR(:,:,1) = eye(global_params.N);
D_AR(:,:,global_params.delay+1) = A_1;
y = customFilter(D_AR,x,global_params.N);
Y = y';
plot(y)

%interval length of 100
interval_length = 50;
intervals = floor(global_params.l/(interval_length/2));
N = global_params.N;
gamma = donoghue_params.gamma;
Y_int = zeros(global_params.N,interval_length,floor(global_params.l/(interval_length/2)));
for i=1:N
    Y_int(i,:,:) = buffer(Y(i,:),interval_length,interval_length/2);
end
L = zeros(N,N,intervals);
for i = 1:intervals
    Y_tmp = Y_int(:,:,i);
    B = approxAR(Y_tmp',donoghue_params.P);
    Q = create_Q_matrix(N,donoghue_params.gamma);
    c = create_c_vec(N,B');
    A = create_constraint_matrix(N);
    b = create_constraint_vec(N,signalEnergy(Y_tmp)*N);
    model.Q = sparse(Q);
    model.A = sparse(A);
    model.obj = c;
    model.rhs = b;
    model.sense = [char('='*ones(N+1,1));char('<'*ones(N*(N-1)/2,1))];
    model.lb = -inf*ones(N*(N+1)/2,1);
    params.outputflag = 0;
    results = gurobi(model,params);
    phi = results.x;
    l = create_dup_matrix(N)*phi;
    L(:,:,i) = convert_to_matrix(l);
end

A_donoghue_occ_t = laplacian_to_adjacency(L)>donoghue_params.threshold; % binarize each interval's A
A_donoghue_occ_avg = mean(A_donoghue_occ_t,3); % average binarized A's
A_donoghue_occ_avg(A_donoghue_occ_avg<donoghue_params.threshold_avg) = 0; % threshold the averaged result
L_donoghue_occ_avg = adjacency_to_laplacian(A_donoghue_occ_avg); % final laplacian


[p,r,f,nm,ne] = graph_learning_perf_eval(L_0,L_donoghue_occ_avg)
vcompare(L_0,L_donoghue_occ_avg)