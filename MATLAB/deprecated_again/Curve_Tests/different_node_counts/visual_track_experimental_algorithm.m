clc;clear;close all hidden;
load('large_maximum_data_results_n100.mat');
load('large_maximum_data_results_n100_newP.mat');
load('large_maximum_data_results_n15_n10_3.mat');
load('large_maximum_data_results_n20_3.mat');
addpath(genpath('../../GL_classes/'));
%%
% N = 20;
% 
% sparsities = S_RND_20;
% sparse_F = F_RND_20;
% graph_type = 'gaussian';
% graph_param = rnd_param;
% 
% sparsities = S_ER_20;
% sparse_F = F_ER_20;
% graph_type = 'er';
% graph_param = er_param;
% 
% sparsities = S_BA_20;
% sparse_F = F_BA_20;
% graph_type = 'pa';
% graph_param = ba_param;

N = 15;

sparsities = S_RND_15;
sparse_F = F_RND_15;
graph_type = 'gaussian';
graph_param = rnd_param;
% 
% sparsities = S_BA_15;
% sparse_F = F_BA_15;
% graph_type = 'pa';
% graph_param = ba_param;
% 
% sparsities = S_ER_15;
% sparse_F = F_ER_15;
% graph_type = 'er';
% graph_param = er_param;


% N = 10;
% 
% sparsities = S_RND_10;
% sparse_F = F_RND_10;
% graph_type = 'gaussian';
% graph_param = rnd_param;

% sparsities = S_BA_10;
% sparse_F = F_BA_10;
% graph_type = 'pa';
% graph_param = ba_param;

% sparsities = S_ER_10;
% sparse_F = F_ER_10;
% graph_type = 'er';
% graph_param = er_param;


% N = 100
% sparsities = S_ER_100;
% sparse_F = F_ER_100;
% graph_type = 'er';
% graph_param = .2;


% thresholds = 0:.0025:2;
thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);

calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(N*(N-1)/2);
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

    % edge_spread = graphs.get_edge_spread(graph_type,graph_param(i),N);
    % edge_spread = graphs.get_edge_spread2(graph_type,graph_param(i),N);
    % v_std = graphs.get_vertice_std(graph_type,graph_param(i),N);
    % sparsity = graphs.get_sparsity(graph_type,graph_param(i));

    % [a,b] = GL.get_cliff_curve(sparsities(i));
    % [c1,c2] = GL.plot_cliff_curve(gammas,a,b,N,v_std);

    % [X,Y] = meshgrid(log10(gammas),thresholds*N);
    % cliff_curve_z1 = interp2(X,Y,sparse_F(:,:,i)',log10(gammas),c1)+1e-3;
    % [X,Y] = meshgrid(log10(gammas),thresholds*N);
    % cliff_curve_z2 = interp2(X,Y,sparse_F(:,:,i)',log10(gammas),c2)+1e-3;

    % plot3(log10(gammas),c1,cliff_curve_z1,'LineWidth',2);
    % plot3(log10(gammas),c2,cliff_curve_z2,'LineWidth',2);
    % plot3(log10(gammas),c2,cliff_curve_z2*10,'LineWidth',2,'Color','red');
    % plot3(log10(gammas),c1,cliff_curve_z1*10,'LineWidth',2);
    % 
    % 
    % t = GL.get_threshold(edge_spread,N)*(N);
    % 
    % plot3(log10(gammas),t*ones(length(gammas),1),ones(length(gammas),1),'LineWidth',2,'Color','black');

    % a = 1.4753;
    % b = -3.4686;
    % c = 1.6991;
    % d = 0.0169;
    % g = a+b*exp(-c*(t-d)).*cos(-c*(t-d));

    % g = GL.get_gamma(sparsities(i),t/N,N,v_std);
    % g = log10(g);
    % 
    % plot3(g*ones(length(thresholds),1),thresholds*N,ones(length(thresholds),1),'LineWidth',2,'Color','black');
    % 
    % [z_all(i),ind] = max(sparse_F(:,:,i),[],"all");
    % [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    % z_g = log10(gammas(row));
    % z_t = thresholds(col)*N;
    % scatter3(z_g,z_t,1);

    % [z, ind] = max(cliff_curve_z1,[],"all");
    % z_g = log10(gammas(ind));
    % z_t = c1(ind);
    % scatter3(z_g,z_t,1,'filled','MarkerFaceColor','black');
    % 
    % rounded_g = interp1(gammas,gammas,10^g,'nearest');
    % rounded_t = interp1(thresholds,thresholds,t/N,'nearest');
    % ind_g = find(gammas == rounded_g);
    % ind_t = find(thresholds == rounded_t);
    % curveF(i) = sparse_F(ind_g,ind_t,i)*100;
    % 
    % 
    % title(sprintf("maxF = %.1d, curveF = %.1d",z_all(i)*100,curveF(i)),sparsity)
    % pause(1)
end
%%
figure;hold on;
plot(z_all*100,'--')
plot(curveF)
ylim([0 100])