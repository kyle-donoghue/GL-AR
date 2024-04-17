clc;clear;close all hidden;
addpath(genpath('../../GL_classes/'));
load('temperature_dataset.mat')
%%
[temperature,mu,sigma] = signals.z_score(temperature); % z-score to help fit
for i = 1:size(temperature,1)
    fittedModel = fit((1:365)'-round(365/2),temperature(i,:)', 'poly4'); % translate to be around 0
    fittedValues(i,:) = feval(fittedModel,(1:365)-round(365/2));
end
fittedValues = fittedValues.*sigma+mu; %inverse z-score to return to normal
temperature = temperature.*sigma+mu;

detrended = temperature-fittedValues;
plot(detrended')
%%
