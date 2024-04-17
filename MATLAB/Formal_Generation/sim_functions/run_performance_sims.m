function [F,densities,T] = run_performance_sims(signal_params,AR_params,graph_params,trials)

    F = zeros(length(graph_params),5);
    T = zeros(length(graph_params),2);
    densities = zeros(length(graph_params),1);
    for s = 1:length(graph_params)
        fprintf("%d %s %d\n",signal_params.N,signal_params.graph_type,graph_params(s))
        toc
        tic
        p = graph_params(s);
        densities(s) = graphs.get_density(signal_params.graph_type,p,signal_params.N);
        density = densities(s);

        AR_params.gamma = GL.get_gamma(density);


        max_edges = graphs.max_edges(signal_params.N);
        num_edges = ceil(density*max_edges);
        pred_edges = max(1,num_edges-1);

        f = zeros(trials,5);
        diff_t = zeros(trials,1);
        parfor i = 1:trials
            [L_0,~,A_d] = graphs.create(signal_params,p);
            G = graphs.createGraphTensor(signal_params,A_d);
            y_noisy = signals.generateFilteredRectPulse(signal_params,G);
            % y_noisy = signals.generateFilteredARProcess(signal_params,G);
            % y_noisy = signals.generateFilteredRectPulseWindow(signal_params,G);
            % y_noisy = signals.generateFilteredSine(signal_params,G);
            % y_noisy = signals.generateFilteredGaussian(signal_params,G);
            % nnz(~y_noisy)
            % nnz(isnan(y_noisy))
            % y_noisy = signals.generateFilteredSinePulse(signal_params,G,2);
            % y_noisy = signals.generateFilteredSinePulseWindow(signal_params,G);
            % y_noisy = signals.generateFilteredNyquistPulse(signal_params,G);
            y_noisy = signals.z_score(y_noisy);
            % nnz(isnan(y_noisy))
            L = GL.AR_mean(y_noisy,AR_params);
    
            weights = graphs.get_weights(L);
            t = GL.get_threshold(weights,pred_edges,density);
            diff_t(i) = t;

            L_tmp = GL.threshold(L,t);
            [f1,f2,f3,f4,f5] = graphs.performance(L_0,L_tmp);
            f(i,:) = [f1 f2 f3 f4 f5];
        end
        F(s,:) = mean(f,1);
        T(s,1) = mean(diff_t);
        T(s,2) = pred_edges;
    end
end
