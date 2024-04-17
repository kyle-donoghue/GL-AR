clc;clear;close all hidden;
%%
load('raw_temperature_data.mat')
ID_array = [];
temperature_signals = NaN(365);
count = 1;
prevInd =- 1;
for i = 1:length(temperature_data)
    ind1 = temperature_data(i,2);
    if ind1 ~= prevInd
        ID_array(end+1) = ind1;
        count = count+1;
        prevInd = ind1;
    end
    ind2 = day(datetime(num2str(temperature_data(i,3)),'InputFormat','yyyyMMdd'),'dayofyear');
    temperature_signals(count,ind2) = temperature_data(i,4);
end
good_temperature_signals = [];
good_ID = [];
count = 1;
for i = 1:size(temperature_signals,1)
    if ~isnan(temperature_signals(i,:))
        good_ID(count) = ID_array(i);
        good_temperature_signals(count,:) = temperature_signals(i,:);
        count = count+1;
    end
end