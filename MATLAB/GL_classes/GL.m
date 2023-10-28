classdef GL
    %GL Class for graph learning algorithms including Dong and AR
    
    methods (Static=true)
        function L = AR(Y,p)
            % bad_data = max(isnan(Y),[],"all");
            
            Q = create_Q_matrix(p.N,1);
            B = signals.approxAR(Y,p.P);
            c = create_c_vec(p.N,B)*p.gamma;
            A = create_constraint_matrix(p.N);
            b = create_constraint_vec(p.N,p.N);

            % options = optimset('Display', 'off');
            % [phi,f] = quadprog(Q,c,A(p.N+2:end,:),b(p.N+2:end),A(1:p.N+1,:),b(1:p.N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
   
            lowerbound = -Inf(size(b));
            lowerbound(1:p.N+1) = b(1:p.N+1);
            upperbound = b;
            m = osqp;
            m.setup(Q,c,A,lowerbound,upperbound,'verbose',false);
            results = m.solve();
            phi = results.x;
            
            l = create_dup_matrix(p.N)*phi;
            L = convert_to_matrix(l);

        end
        function [L,Y,L_harvard] = dong(X_noisy,param)
            N = param.N;
            max_iter = param.max_iter;
            alpha = param.alpha;
            beta = param.beta;
            objective = zeros(max_iter,1);
            Y_0 = X_noisy;
            Y = Y_0;
            for i = 1:max_iter
                % Step 1: given Y, update L
                L = optimize_laplacian_gaussian(N,Y,alpha,beta);
            %     L = optimize_laplacian_gaussian_admm(N,Y,alpha,beta,0.1,1.5);
                % solution in the Harvard paper
                if i == 1
                    L_harvard = L;
                end
                % Step 2: Given L, update Y
                % Y = (eye(N)+alpha*L)^(-1)*Y_0;
                R = chol(eye(N) + alpha*L);
                Y = R \ (R' \ (Y_0));
                % plot the objective
                % objective(i) = norm(Y-Y_0,'fro')^2 + alpha*trace(Y'*L*Y) + beta*(norm(L,'fro')^2);
                objective(i) = norm(Y-Y_0,'fro')^2 + alpha*vec(Y*Y')'*vec(L) + beta*(norm(L,'fro')^2);
                % stopping criteria
                if i>=2 && abs(objective(i)-objective(i-1))<10^(-4)
                    break
                end
            end
        end
        function L = threshold(L,t)
            L(abs(L)<=t) = 0;
        end
    end
end

