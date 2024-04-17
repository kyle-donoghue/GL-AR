clc;clear;close all hidden;
load('thresholds_n20.mat');
load('thresholds_n10.mat');
% load('large_maximum_data_results_n20_3.mat');
% load('large_maximum_data_results_n15_n10_3.mat');
% load('large_maximum_data_results_count_er_100_increasedSNR2.mat');
% load('large_maximum_data_results_count_er_signaltest_SNR200.mat');
% load('large_maximum_data_results_count_er_ztest2.mat');
% load('large_maximum_data_results_count_er_ztest200.mat');
% load('large_maximum_data_results_count_er_ztest200_N20.mat');

% load("large_new_performance_data_results_z_n40.mat")
load('large_new_maximum_data_results_n40_n20_n10.mat');
addpath(genpath('../../GL_classes/'));
%%

% N = 20;
% sparsities = S_RND_20;
% sparse_F = F_RND_20;
% graph_type = 'gaussian';
% graph_param = rnd_param;
% T = T_RND_20;

% N = 20;
% sparsities = S_BA_20;
% sparse_F = F_BA_20;
% graph_type = 'pa';
% graph_param = ba_param;
% T = T_BA_20;

N = 20;
sparsities = S_ER_20;
sparse_F = F_ER_20;
graph_type = 'er';
graph_param = er_param;
T = T_ER_20;

% N=10;
% sparsities = S_ER_10;
% sparse_F = F_ER_10;
% graph_type = 'er';
% graph_param = er_param;
% T = T_ER_10;

% N=10;
% sparsities = S_RND_10;
% sparse_F = F_RND_10;
% graph_type = 'gaussian';
% graph_param = rnd_param;
% T = T_RND_10;

% N=10;
% sparsities = S_BA_10;
% sparse_F = F_BA_10;
% graph_type = 'pa';
% graph_param = ba_param;
% T = T_BA_10;\

% N=40;
% sparsities = S_ER_40;
% sparse_F = F_ER_40;
% graph_type = 'er';
% graph_param = .1:.1:.9;
% T = T_ER_40;

% N=N_list(1);
% sparsities = .2;
% sparse_F = F_surfs(:,:,1);
% graph_type = 'er';
% graph_param = .2;
% T = T_ER_20;

thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
for i = 10:(length(sparsities))
    figure(1);
    clf;hold on;
    fig.fig = gcf;
    [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds(2:end)*N);
    fig.s = surf(X_axis,Y_axis,sparse_F(:,2:end,i)');
    fig.s.EdgeColor = "none";
    xlabel("\textbf{log}$$_{\textbf{10}}(\gamma)$$",'Interpreter','latex','FontSize',14)
    ylabel("$$\tilde{\delta}$$",'Interpreter','latex','FontSize',14)
    zlabel("\textbf{F-score}",'Interpreter','latex','FontSize',14)
    view(2)
    ylim([0,10])    
    STANDARDIZE_FIGURE(fig)

    % edge_spread = graphs.get_edge_spread(graph_type,graph_param(i),N);
    edge_spread = graphs.get_edge_spread2(graph_type,graph_param(i),N);
    v_std = graphs.get_vertice_std(graph_type,graph_param(i),N);
    sparsity = graphs.get_sparsity(graph_type,graph_param(i));

    [a,b] = GL.get_cliff_curve(sparsities(i));
    [c1,c2] = GL.plot_cliff_curve(gammas,a,b,N,v_std);

    [X,Y] = meshgrid(log10(gammas),thresholds*N);
    cliff_curve_z1 = interp2(X,Y,sparse_F(:,:,i)',log10(gammas),c1)+1e-3;
    [X,Y] = meshgrid(log10(gammas),thresholds*N);
    cliff_curve_z2 = interp2(X,Y,sparse_F(:,:,i)',log10(gammas),c2)+1e-3;

    plot3(log10(gammas),c1,cliff_curve_z1,'LineWidth',2,'Color','red');
    plot3(log10(gammas),c2,cliff_curve_z2,'LineWidth',2);

    % a = -67.8256;
    % b = 0.0413;
    % c = 69.3318;
    % g = c+a*exp(b.*sparsity.^2);
    a = -0.2724;
    b = 2.6540;
    g = a*exp(b.*sparsity);
    rounded_g = interp1(gammas,gammas,10^g,'nearest');
    ind_g = find(gammas == rounded_g);
    % plot3(g*ones(length(thresholds),1),thresholds*N,max(sparse_F(ind_g,:,i),[],"all")*ones(length(thresholds),1),'LineWidth',2,'Color','black');
    plot3(g*ones(length(thresholds),1),thresholds*N,max(sparse_F(:,:,i),[],"all")*ones(length(thresholds),1),'LineWidth',2,'Color','black');
    
    t = T(i,1)*N;
    %make ramp function from (.3,0) to (.7,.25)
    % if sparsity >= .25 && sparsity <= .55
    %     ramp = (sparsity-.25)
    %     t = t-ramp;
    % elseif sparsity > .55
    %     ramp = .3
    %     t = t-ramp;
    % else
    %     ramp = 0
    % end

    plot3(log10(gammas),t*ones(length(gammas),1),max(sparse_F(:,:,i),[],"all")*ones(length(gammas),1),'LineWidth',2,'Color','black');

    [z_all(i),ind] = max(sparse_F(:,:,i),[],"all");

    % [z, ind] = max(cliff_curve_z1,[],"all");
    % z_g = log10(gammas(ind));
    % z_t = c1(ind);
    % scatter3(z_g,z_t,1,'filled','MarkerFaceColor','black');

    rounded_g = interp1(gammas,gammas,10^g,'nearest');
    rounded_t = interp1(thresholds,thresholds,t/N,'nearest');
    ind_g = find(gammas == rounded_g);
    ind_t = find(thresholds == rounded_t);
    curveF(i) = sparse_F(ind_g,ind_t,i)*100;


    title(sprintf("maxF = %.1d, curveF = %.1d",z_all(i)*100,curveF(i)),sparsity)


    % pause(1)
end