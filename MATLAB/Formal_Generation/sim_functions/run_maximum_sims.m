function [sparse_F,densities] = run_maximum_sims(signal_params,graph_params,trials,gammas,thresholds)
    AR_params = GL.create_default_params(signal_params);
    
    t_max = length(thresholds);

    sparse_F = zeros(length(gammas),length(thresholds),length(graph_params));
    densities = zeros(length(graph_params),1);
    count = 1;
    for s = 1:length(graph_params)
        p = graph_params(s);
        densities(s) = graphs.get_density(signal_params.graph_type,p);
        F = zeros(length(gammas),length(thresholds));
        for k = 1:length(gammas)
            fprintf("%d %s",signal_params.N,signal_params.graph_type)
            tim = toc
            tic
            percent = count/(length(gammas)*length(graph_params))*100
            left = tim*((length(gammas)*length(graph_params))-count)/60
            count = count+1;
    
            
            AR_params.gamma = gammas(k);
        
            f = zeros(t_max,trials);
            parfor i = 1:trials
                [L_0,~,A_d] = graphs.create(signal_params,p);
                G = graphs.createGraphTensor(signal_params,A_d);
                % y_noisy = signals.generateFilteredRectPulse(signal_params,G);
                % y_noisy = signals.generateFilteredARProcess(signal_params,G);
                y_noisy = signals.generateFilteredSine(signal_params,G);
                y_noisy = signals.z_score(y_noisy);

                L = GL.AR_mean(y_noisy,AR_params);
        
                for t = 1:t_max 
                    L_tmp = GL.threshold(L,thresholds(t));
                    [~,~,f(t,i),~,~] = graphs.performance(L_0,L_tmp);
                end
            end
            F(k,:) = mean(f,2);
        end
        sparse_F(:,:,s) = F;
        max(F,[],"all")
    end
    max_F = max(sparse_F,[],"all");
end