clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));

%%
dataset = readtable('ad_viz_plotval_data.csv');
dataset.Date = day(dataset.Date,'dayofyear');
index = 0;
prevID = -1;
for i = 1:height(dataset)
    currID = dataset.SiteID(i);
    if currID ~= prevID
        index = index+1;
    end
    prevID = currID;
    day = dataset.Date(i);
    pollution(index,day) = dataset.DailyMeanPM2_5Concentration(i);
end
poll_energy = sum(pollution.^2,2);
[poll_energy,ind] = sort(poll_energy,'descend');
plot(pollution(ind(1:20),:)')
pollution = pollution(ind(1:20),:);

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
latlong = latlong(ind(1:20),:);
A = zeros(20);
for i = 1:20
    for j = (i+1):20
        A(i,j) = lldistkm(latlong(i,:),latlong(j,:));
    end
end
dist_adj = A+A';
figure;
heatmap(dist_adj)
save("pollution_data.mat",'pollution','dist_adj')