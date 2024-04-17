clc;clear;close all hidden;
addpath(genpath("../comparison_Tests/"));
rng("default")
global_params.N = 20; % node count
global_params.l = 1e3; % signal length
global_params.noise = .000005; %noise power
global_params.trials = 50;
global_params.Fs = 100;
global_params.delay = 10;
global_params.minSep = 100;
global_params.maxSep = 250;
donoghue_params = global_params;
donoghue_params.P = 15; % AR approx order
donoghue_params.gamma = .002;
donoghue_params.eig_padding = .99; % what should the max eigenvalue of A be
donoghue_params.intervals = 20;
donoghue_params.interval_length = floor(global_params.l/donoghue_params.intervals);
donoghue_params.threshold_avg = 0;
donoghue_params.threshold = 10.5;

[A,XCoords, YCoords] = construct_graph(global_params.N,'er',0.2);
A_0 = A*(donoghue_params.eig_padding/max(abs(eig(A)))); % make VARM stable by normalizing eigenvalues
L_0 = adjacency_to_laplacian(A_0);
A_1 = to_directed(A_0);
freqs = 10.2:.2:14;
freqs = generate_frequencies2(global_params.N,2,25);

x = zeros(global_params.l,global_params.N);
for i = [1:global_params.N]
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

%interval length of 100
interval_length = 50;
intervals = floor(global_params.l/(interval_length/2));
N = global_params.N;
Y_int = zeros(global_params.N,interval_length,floor(global_params.l/(interval_length/2)));
for i=1:N
    Y_int(i,:,:) = buffer(Y(i,:),interval_length,interval_length/2);
end
figure
plot(Y')
figure
imagesc(A_1)
Y = Y_int(:,:,6);
Y2 = Y_int(:,:,11);
Y3 = Y_int(:,:,20);





[e1, e1_i] = signalEnergy(Y);
[e2, e2_i] = signalEnergy(Y2);
[e3, e3_i] = signalEnergy(Y3);










B = approxAR(Y',donoghue_params.P);
Q = create_Q_matrix(N,donoghue_params.gamma);
c = create_c_vec(N,B(2:end,:)');
A = create_constraint_matrix(N);
b = create_constraint_vec(N,sqrt(e1)*N);
b(N+1) = e1*N;
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
L = convert_to_matrix(l);
A = laplacian_to_adjacency(L);
A(A<donoghue_params.threshold) = 0;
L = adjacency_to_laplacian(A);

B2 = approxAR(Y2',donoghue_params.P);
Q = create_Q_matrix(N,donoghue_params.gamma);
c = create_c_vec(N,B2(2:end,:)');
A = create_constraint_matrix(N);
b = create_constraint_vec(N,sqrt(e2)*N);
b(N+1) = e2*N;
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
L2 = convert_to_matrix(l);
A = laplacian_to_adjacency(L2);
A(A<donoghue_params.threshold) = 0;
L2 = adjacency_to_laplacian(A);

B3 = approxAR(Y3',donoghue_params.P);
Q = create_Q_matrix(N,donoghue_params.gamma);
c = create_c_vec(N,B3(2:end,:)');
A = create_constraint_matrix(N);
b = create_constraint_vec(N,sqrt(e3)*N);
b(N+1) = e3*N;
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
L3 = convert_to_matrix(l);
A = laplacian_to_adjacency(L3);
A(A<donoghue_params.threshold) = 0;
L3 = adjacency_to_laplacian(A);

[p,r,f,nm,ne] = graph_learning_perf_eval(L_0,L);

figure;
subplot(1,3,1)
plot(Y')
legendStrings = "N = " + string(1:N);
legend(legendStrings)
subplot(1,3,2)
plot(Y2')
legendStrings = "N = " + string(1:N);
legend(legendStrings)
subplot(1,3,3)
plot(Y3')
legendStrings = "N = " + string(1:N);
legend(legendStrings)

figure
subplot(1,3,1)
plot(B(2:end,:))
legendStrings = "N = " + string(1:N);
legend(legendStrings)
subplot(1,3,2)
plot(B2(2:end,:))
legendStrings = "N = " + string(1:N);
legend(legendStrings)
subplot(1,3,3)
plot(B3(2:end,:))
legendStrings = "N = " + string(1:N);
legend(legendStrings)

figure
subplot(1,2,1)
plot(abs(freqz(1,B(:,1))));
subplot(1,2,2)
plot(abs(freqz(1,B2(:,N))));
vcompare(L_0,L);
vcompare(L_0,L2);
vcompare(L_0,L3);


