clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('temperature_dataset.mat')
%%
good_indices = [125 24 25 126 23 124 115 70 130 7 87 49 119 138 21 97 20 93 57 145 63 94 14 16  88 103 73 82 39 72 141 96 70 17 149 117 59 145 52 130 86 135 78 64 34 35 104 67 102 46];
latlong = cell2mat(cities(:,2:3));
legend_list = {};
figure;hold on;
for i = 1:length(DIST)
    scatter(latlong(i,1),latlong(i,2))
    legend_list{end+1} = num2str(i);
end
legend(legend_list)