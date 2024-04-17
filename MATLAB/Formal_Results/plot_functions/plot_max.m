function p = plot_max(i,j,s,f,M)
    figure(i);hold on;
    % subplot(1,M,j);hold on;
    p = plot(s,permute(max(f,[],[1 2]),[3 1 2]),'LineWidth',3);
    xlabel('\textbf{Graph Size}','Interpreter','latex','FontSize',14)
    ylabel('\textbf{F-score}','Interpreter','latex','FontSize',14)
    % title('ER','Interpreter','latex','FontSize',14)
    % legend('Algorithm Performance','Theoretical Performance')
    % ylim([0 109])
end