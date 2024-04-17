clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
%%
load('large_maximum_data_results_count_er_100_increasedSNR2.mat')

index = 14; %N=46
N = N_list(index);
sparsity = .2;
graph_param = .2;
graph_type = 'er';



thresholds = 0:.0025:2;
gammas = logspace(-2,4,60);
%%
figure(1);
clf;hold on;
[X_axis,Y_axis] = meshgrid(log10(gammas),thresholds*N);
s = surf(X_axis,Y_axis,F_surfs(:,:,index)');
s.EdgeColor = "none";
xlabel("log$$_{10}(\gamma)$$",'Interpreter','latex','FontSize',14)
ylabel("$$\tilde{\delta}$$",'Interpreter','latex','FontSize',14)
zlabel("F-score",'Interpreter','latex','FontSize',14)
view(2)
ylim([0,20])    

% edge_spread = graphs.get_edge_spread(graph_type,graph_param(i),N);
edge_spread = graphs.get_edge_spread2(graph_type,graph_param,N);
v_std = graphs.get_vertice_std(graph_type,graph_param,N);
sparsity = graphs.get_sparsity(graph_type,graph_param);

[a,b] = GL.get_cliff_curve(sparsity);
[c1,c2] = GL.plot_cliff_curve(gammas,a,b,N,v_std);

[X,Y] = meshgrid(log10(gammas),thresholds*N);
cliff_curve_z1 = interp2(X,Y,F_surfs(:,:,index)',log10(gammas),c1)+1e-3;
[X,Y] = meshgrid(log10(gammas),thresholds*N);
cliff_curve_z2 = interp2(X,Y,F_surfs(:,:,index)',log10(gammas),c2)+1e-3;

plot3(log10(gammas),c1,cliff_curve_z1,'LineWidth',2,'Color','red');
plot3(log10(gammas),c2,cliff_curve_z2,'LineWidth',2);


t = GL.get_threshold(edge_spread,N)*(N);
g = GL.get_gamma(sparsity,t/N,N,v_std);
g = log10(g);

plot3(log10(gammas),t*ones(length(gammas),1),ones(length(gammas),1),'LineWidth',2,'Color','black');
plot3(g*ones(length(thresholds),1),thresholds*N,ones(length(thresholds),1),'LineWidth',2,'Color','black');

[z_all,~] = max(F_surfs(:,:,index),[],"all");

[z, ind] = max(cliff_curve_z1,[],"all");
z_g = log10(gammas(ind));
z_t = c1(ind);
scatter3(z_g,z_t,1,'filled','MarkerFaceColor','black');

rounded_g = interp1(gammas,gammas,10^g,'nearest');
rounded_t = interp1(thresholds,thresholds,t/N,'nearest');
ind_g = find(gammas == rounded_g);
ind_t = find(thresholds == rounded_t);
curveF = F_surfs(ind_g,ind_t,index)*100;


title(sprintf("maxF = %.1d, curveF = %.1d",z_all*100,curveF),sparsity)