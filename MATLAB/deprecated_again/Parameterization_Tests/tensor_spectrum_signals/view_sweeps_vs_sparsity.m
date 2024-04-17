clc;clear;close all hidden;
load("gamma_threshold_sweep_with_sparsity_small_BA.mat");
%%
max_gammas = zeros(length(sparsities),1);
for i = 1:length(sparsities)
    [~, ind] = max(sparse_F(:,:,i),[],"all");
    [row, col] = ind2sub(size(F),ind);
    max_gammas(i) = gammas(row);
end
%%
a_a = 0.1289;
a_b = -5.855;
a_offset = 0.0060;
b_a1 =  0.4011;
b_c1 = 0.2049;
b_offset = 1.5;
cliff_curve_offset = 1/signal_params.N;
cliff_curve_delay =1;
% ridge_offest = .75;
ridge_offset = .1;
%%
% t_a = 0.4285;
% t_b = -12.86;
% t_offset = .75/signal_params.N;
t_a = 0.4285; % changed to be with num_edges instead of sparsity
% t_b = -0.6431;
t_b = -0.06768;
t_offset = .75/signal_params.N;

%%
figure;
calculated = zeros(length(sparsities),1);
maximums = zeros(length(sparsities),1);
num_edges = sparsities*(signal_params.N*(signal_params.N-1)/2);
for i = 1:length(sparsities)
    clf;hold on;
    [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds);
    s = surf(X_axis,Y_axis,sparse_F(:,:,i)');
    s.EdgeColor = "none";
    view(2)

    a = a_a*exp(a_b*sparsities(i))+a_offset;
    b = b_a1*exp(-((sparsities(i))/b_c1)^2)+b_offset;

    cliff_curve = a*exp(b*log10(gammas))+cliff_curve_offset;
    cliff_curve = cliff_curve(cliff_curve <= thresholds(end));
    cliff_curve2 = a*exp(b*(log10(gammas)-cliff_curve_delay))+ridge_offset*cliff_curve_offset;
    cliff_curve2 = cliff_curve2(cliff_curve2 <= thresholds(end));
    plot3(log10(gammas(1:length(cliff_curve))),cliff_curve,ones(length(cliff_curve),1),'LineWidth',2);
    plot3(log10(gammas(1:length(cliff_curve2))),cliff_curve2,ones(length(cliff_curve2),1),'LineWidth',2);
    



    t = t_a*exp(t_b*(num_edges(i)-num_edges(1)))+t_offset;
    
    % if t > cliff_curve_offset/2
    if t > t_offset*1.05
        g = max((log((t-ridge_offset*cliff_curve_offset)/a)/b+cliff_curve_delay),0);
    else
        g = 0;
    end

    rounded_g = interp1(gammas,gammas,10^g,'nearest');
    rounded_t = interp1(thresholds,thresholds,t,'nearest');
    ind_g = find(gammas == rounded_g);
    ind_t = find(thresholds == rounded_t);

    scatter3(g,t,1,'MarkerFaceColor','red')
    scatter3(log10(rounded_g),rounded_t,1,'MarkerFaceColor','magenta')

    [maxF,ind] = max(sparse_F(:,:,i),[],"all");
    curveF = sparse_F(ind_g,ind_t,i)*100;
    [row, col] = ind2sub(size(sparse_F(:,:,i)),ind);
    max_g = gammas(row);
    max_t = thresholds(col);
    title(sprintf("s = %.3d, maxF = %.1d, curveF = %.1d",sparsities(i),maxF*100,curveF))

    scatter3(log10(max_g),max_t,1,'MarkerFaceColor','cyan')

    maximums(i) = maxF*100;
    calculated(i) = curveF;
    xlabel("Gamma")
    ylabel("Threshold")
    % pause(1)
end
%%
% load("calculated_curves.mat")
figure;hold on;
plot(sparsities,maximums,'--');
plot(sparsities,calculated)
ylim([0 100])
xlabel('Sparsity')
ylabel('F-Score')
max(maximums-calculated)
mse(maximums,calculated)