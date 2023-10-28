function tcompare(L1,L2)
    thresh = .01:.01:.3;
    f = zeros(1,length(thresh));
    for t = 1:length(thresh)
        L2tmp = L2;
        L2tmp(abs(L2tmp)<thresh(t)) = 0;
        [~,~,f(t),~,~] = graph_learning_perf_eval(L1,L2tmp);
    end
    figure;
    plot(thresh,f);
end