function [sparse_F,densities,sparse_p,sparse_r,sparse_nmi] = run_maximum_sims_dong(signal_params,graph_params,trials,gammas,thresholds)
    dong_params = signal_params;
    dong_params.l = dong_params.interval_length;
    dong_params.max_iter = 50;
    dong_params.alpha = 10.^(-2);
    dong_params.beta = 10.^(-0.2);
    dong_params.lambda = 10.^1;

    t_max = length(thresholds);

    sparse_F = zeros(length(gammas),length(thresholds),length(graph_params));
    densities = zeros(length(graph_params),1);
    count = 1;
    for s = 1:length(graph_params)
        p = graph_params(s);
        densities(s) = graphs.get_density(signal_params.graph_type,p,signal_params.N);
        F = zeros(length(gammas),length(thresholds));
        PRE = zeros(length(gammas),length(thresholds));
        REC = zeros(length(gammas),length(thresholds));
        NMI = zeros(length(gammas),length(thresholds));
        for k = 1:length(gammas)
            fprintf("Dong %d %s",signal_params.N,signal_params.graph_type)
            tim = toc
            tic
            percent = count/(length(gammas)*length(graph_params))*100
            left = tim*((length(gammas)*length(graph_params))-count)/60
            count = count+1;
    
            
            dong_params.beta = gammas(k);
        
            f = zeros(t_max,trials);
            pre = zeros(t_max,trials);
            rec = zeros(t_max,trials);
            nmi = zeros(t_max,trials);
            parfor i = 1:trials
                [L_0,~,A_d] = graphs.create(signal_params,p);
                G = graphs.createGraphTensor(signal_params,A_d);
                % y_noisy = signals.generateFilteredSine(signal_params,G);
                y_noisy = signals.generateFilteredRectPulse(signal_params,G);
                % y_noisy = signals.generateFilteredSinePulseWindow(signal_params,G);
                y_noisy = signals.z_score(y_noisy);

                [L,~,~] = GL.dong(y_noisy(:,1:signal_params.interval_length),dong_params);
        
                for t = 1:t_max 
                    L_tmp = GL.threshold(L,thresholds(t));
                    [pre(t,i),rec(t,i),f(t,i),nmi(t,i),~] = graphs.performance(L_0,L_tmp);
                end
            end
            F(k,:) = mean(f,2);
            PRE(k,:) = mean(pre,2);
            REC(k,:) = mean(rec,2);
            NMI(k,:) = mean(nmi,2);
        end
        sparse_F(:,:,s) = F;
        sparse_p(:,:,s) = PRE;
        sparse_r(:,:,s) = REC;
        sparse_nmi(:,:,s) = NMI;
    end
    max_F = max(sparse_F,[],"all");
end