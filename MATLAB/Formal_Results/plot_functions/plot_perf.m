function p = plot_perf(i,j,s,f,M)
    figure(i);hold on;
    subplot(1,M,j);hold on;
    p = plot(s,f,'--','LineWidth',3);
    ylim([0 1.1])
    % xlim([.5 10.5])
end