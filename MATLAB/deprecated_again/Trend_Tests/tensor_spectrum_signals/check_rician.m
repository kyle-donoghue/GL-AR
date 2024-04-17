clc;clear;close all hidden;
% A = [.1:.05:2-.05];
% B = .5*2.^(1:.05:2 -.05);
% B = [B./(1.95:-.05:1) fliplr(B)./(1:.05:1.95)];
% x = 0:.001:10;
% figure;
% for i = 1:length(A)
%     y = pdf('Rician',x,A(i),B(i));
%     plot(x,y);
%     title('Gamma is decreasing (Sparsity is decreasing) (NumEdges is Increasing)')
%     ylim([0,1])
%     xline(x*y'/1000)
%     pause(.2)
% end

A = [.1:.05:2-.05];
x = 0:.001:10;
figure;
for i = 1:length(A)
    sigma = findSigmaForMean(sqrt(A(i)), 2, 5);
    y = pdf('Rician',x,A(i),sigma^2);
    plot(x,y);
    title('Gamma is decreasing (Sparsity is decreasing) (NumEdges is Increasing)')
    ylim([0,1])
    xline(x*y'/1000)
    pause(.2)
end

function sigma = findSigmaForMean(v, targetMean, initialGuessSigmaSquared)
    % Define the function for the equation to be solved
    meanEquation = @(sigmaSquared) sigmaSquared * sqrt(pi/2) * (1 + v^2 / (2 * sigmaSquared)) * exp(-v^2 / (2 * sigmaSquared)) - targetMean;

    % Set a default initial guess if not provided
    if nargin < 3
        initialGuessSigmaSquared = 1;
    end

    % Solve for sigmaSquared
    sigmaSquared = fsolve(meanEquation, initialGuessSigmaSquared);

    % Calculate sigma from sigmaSquared
    sigma = sqrt(sigmaSquared);
end
