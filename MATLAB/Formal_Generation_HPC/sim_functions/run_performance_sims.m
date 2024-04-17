function [F,densities,T] = run_performance_sims(signal_params,AR_params,graph_params,trials)

    F = zeros(length(graph_params),1);
    T = zeros(length(graph_params),2);
    densities = zeros(length(graph_params),1);
    for s = 1:length(graph_params)
        fprintf("%d %s",signal_params.N,signal_params.graph_type)
        toc
        tic
        p = graph_params(s);
        densities(s) = graphs.get_density(signal_params.graph_type,p,signal_params.N);
        density = densities(s);

        AR_params.gamma = GL.get_gamma(density);


        max_edges = graphs.max_edges(signal_params.N);
        num_edges = ceil(density*max_edges);
        pred_edges = max(1,num_edges-1);

        f = zeros(trials,1);
        diff_t = zeros(trials,1);
        parfor i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,p);
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
            y_noisy = signals.z_score(y_noisy);

            L = GL.AR_mean(y_noisy,AR_params);
    
            weights = graphs.get_weights(L);
            t = GL.get_threshold(weights,pred_edges);
            diff_t(i) = t;

            L_tmp = GL.threshold(L,t);
            [~,~,f(i),~,~] = graphs.performance(L_0,L_tmp);
        end
        F(s) = mean(f);
        T(s,1) = mean(diff_t);
        T(s,2) = pred_edges;
    end
end
