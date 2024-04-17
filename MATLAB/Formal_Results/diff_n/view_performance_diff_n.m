clc;clear;close all hidden;
addpath('../plot_functions');
legend_list = {};
%%
% load('../../Formal_Data/diff_n/performance_data_diff_n_ER_sinepulsewindow4.mat')
% plot_perf(1,1,N_list,F_list(:,3),2)
% defaults(signal_params.graph_type,'Fscore');
% plot_perf(1,2,N_list,F_list(:,4),2)
% defaults(signal_params.graph_type,'NMI');
%%
% load('../../Formal_Data/diff_n/performance_data_diff_n_ER_sine.mat')
% plot_perf(1,1,N_list,F_list(:,3),2)
% defaults(signal_params.graph_type,'Fscore');
% plot_perf(1,2,N_list,F_list(:,4),2)
% defaults(signal_params.graph_type,'NMI');
%%
figure;
fig.fig = gcf;
%%
load('../../Formal_Data/diff_n/performance_data_diff_n_RND_sine_FINAL.mat')
fig.p1 = plot_perf(1,1,N_list,F_list(:,3),1)
defaults(signal_params.graph_type,'Fscore');
% plot_perf(2,2,N_list,F_list(:,4),2)
% defaults(signal_params.graph_type,'NMI');
%%
legend_list{end+1} = 'GL-AR, E$[d]=125$';
%%
load('../../Formal_Data/diff_n/maximum_data_diff_n_RND_dong_sine_FINAL_nodelay.mat')
fig.p2 = plot_max(1,1,N_list,F_surfs,1);
legend_list{end+1} = 'GL-SigRep, E$[d]=0$';
 load('../../Formal_Data/diff_n/maximum_data_diff_n_RND_dong_sine_FINAL_1delay.mat')
fig.p3 = plot_max(1,1,N_list,F_surfs,1);
legend_list{end+1} = 'GL-SigRep, E$[d]=1$';
load('../../Formal_Data/diff_n/maximum_data_diff_n_RND_dong_sine_FINAL_2delay.mat')
fig.p4 = plot_max(1,1,N_list,F_surfs,1);
legend_list{end+1} = 'GL-SigRep, E$[d]=2$';
load('../../Formal_Data/diff_n/maximum_data_diff_n_RND_dong_sine_FINAL_5delay.mat')
fig.p5 = plot_max(1,1,N_list,F_surfs,1);
legend_list{end+1} = 'GL-SigRep, E$[d]=5$';
load('../../Formal_Data/diff_n/maximum_data_diff_n_RND_dong_sine_FINAL_10delay.mat')
fig.p6 = plot_max(1,1,N_list,F_surfs,1);
legend_list{end+1} = 'GL-SigRep, E$[d]=10$';
plot(N_list, .35*ones(size(N_list)),'--','color','black')
%%
% plot(N_list,graphs.minFscore(.2)*100*ones(size(N_list)),'--','Color','black')
%%
xlabel('\textbf{N}')
ylabel('\textbf{F-score}')
% title(signal_params.graph_type)
ylim([0 1])
grid on;
l = legend(legend_list, 'Interpreter', 'Latex');
STANDARDIZE_FIGURE(fig)

%%
% load('../../Formal_Data/diff_n/maxkldsjghkadsj')
% plot_max(1,1,N_list,F_surfs,1)
% legend_list{end+1} = 'GL-AR Bruteforce';

function defaults(t,y)
    ylabel(y)
    xlabel('N')
    % title(t)
    % ylim([0 100])
    grid on;
end