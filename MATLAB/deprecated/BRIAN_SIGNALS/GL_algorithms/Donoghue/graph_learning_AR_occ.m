function L = graph_learning_AR_occ(Y, p)
    debug = 0; %% set to 0 if not plotting


    N = p.N;
    gamma = p.gamma;
    intervals = ceil(p.l/(p.interval_length/2));

    Y_int = zeros(N,p.interval_length,intervals);
    for i=1:N
        Y_int(i,:,:) = buffer(Y(i,:),p.interval_length,(p.interval_length/2));
    end

    L = zeros(p.N,p.N,intervals);


    A = create_constraint_matrix(N);
    % A = create_constraint_matrix2(N);
    Q = create_Q_matrix(N,1);



    % Q = create_Q_matrix(N,gamma);
    % model.Q = sparse(Q);
    % model.A = sparse(A);
    % model.sense = [char('='*ones(N+1,1));char('<'*ones(N*(N-1)/2,1))];
    % model.lb = -inf*ones(N*(N+1)/2,1);
    % params.outputflag = 0;

    for i = 1:intervals
        [e,e_i] = signalEnergy(Y_int(:,:,i));
        B = approxAR(Y_int(:,:,i)',p.P);
        c = gamma*create_c_vec(N,B');
        b = create_constraint_vec(N,sqrt(e)*N);
        % b = create_constraint_vec2(N,e_i,sqrt(e)*N);
        
        % model.rhs = b;
        % model.obj = c;
        % 
        % results = gurobi(model,params);
        % try
        %     phi = results.x;
        % catch
        %     figure
        %     subplot(1,2,1);
        %     plot(Y_int(:,:,i)')
        %     subplot(1,2,2);
        %     plot(B);
        %     error = Y_int(:,:,i);
        %     save("error_data.mat","error")
        % end


        % 
        % tic
        % cvx_begin quiet
        % variable l_0(size(Q,2))
        % minimize (1/2*l_0'*Q*l_0+c'*l_0)
        % Aineq = A(N+2:end,:);
        % Aeq = A(1:N+1,:);
        % bineq = b(N+2:end);
        % beq = b(1:N+1);
        % subject to 
        % Aeq*l_0 == beq
        % Aineq*l_0 <= bineq
        % cvx_end
        % toc

        options = optimset('Display', 'off');
        
        [phi,f] = quadprog(Q,c,A(N+2:end,:),b(N+2:end),A(1:N+1,:),b(1:N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
        % [phi,f] = quadprog(Q,c,A(N+1:end,:),b(N+1:end),A(1:N,:),b(1:N),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
        
        l = create_dup_matrix(N)*phi;
        L(:,:,i) = convert_to_matrix(l);
        if debug
            subplot(2,2,1)
            signalPlot(Y_int(:,:,i)')
            subplot(2,2,2)
            signalPlot(B)
            subplot(2,2,3)
            L_tmp = L(:,:,i);
            matrixPlot(L_tmp)
            subplot(2,2,4)
            L_tmp(abs(L_tmp) < p.threshold) = 0;
            matrixPlot(L_tmp)
        end
    end
end