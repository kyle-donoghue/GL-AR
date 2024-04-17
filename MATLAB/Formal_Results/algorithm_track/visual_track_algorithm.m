clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));

% load('../../Formal_Data/diff_density/maximum_data_diff_density_PARAMETERS.mat');
% load('../../Formal_Data/diff_density/maximum_data_diff_density_default.mat');
load('../../Formal_Data/diff_density/maximum_data_diff_density_L1_default2.mat');
% load('../../Formal_Data/diff_density/performance_data_diff_density_PARAMETERS.mat');
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
% T = T_ER_20;

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


for i = 1:(length(sparsities))
    figure(1);
    clf;hold on;
    [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds*N);
    s = surf(X_axis,Y_axis,sparse_F(:,:,i)');
    s.EdgeColor = "none";
    xlabel("log$$_{10}(\gamma)$$",'Interpreter','latex','FontSize',14)
    ylabel("$$\tilde{\delta}$$",'Interpreter','latex','FontSize',14)
    zlabel("F-score",'Interpreter','latex','FontSize',14)
    view(2)
    ylim([0,20])    

    density = graphs.get_density(graph_type,graph_param(i));

    g = GL.get_gamma(density);

    rounded_g = interp1(gammas,gammas,10^g,'nearest');

    ind_g = (gammas == rounded_g);
    % plot3(g*ones(length(thresholds),1),thresholds*N,max(sparse_F(:,:,i),[],"all")*ones(length(thresholds),1),'LineWidth',2,'Color','black');
    
    % t = T(i,1)*N;

    % plot3(log10(gammas),t*ones(length(gammas),1),max(sparse_F(:,:,i),[],"all")*ones(length(gammas),1),'LineWidth',2,'Color','black');

    [z_all,~] = max(sparse_F(:,:,i),[],"all");
    [z_gamma,~] = max(sparse_F(ind_g,:,i),[],"all");


    % title(sprintf("maxF = %.1d, curveF = %.1d",z_all*100,curveF*100),density)

    % pause(1)
end