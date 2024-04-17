clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));

%%
dataset = readtable('ad_viz_plotval_data.csv');
index = 0;
prevID = -1;
for i = 1:height(dataset)
    currID = dataset.SiteID(i);
    if currID ~= prevID
        index = index+1;
        latlong(index,:) = [dataset.SITE_LATITUDE(i),dataset.SITE_LONGITUDE(i)];
    end
    prevID = currID;
end

for i = 1:20