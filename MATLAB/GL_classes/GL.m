classdef GL
    %GL Class for graph learning algorithms including Dong and AR
    
    methods (Static=true)
        function L = AR(Y,p)            
            Q = create_Q_matrix(p.N,1);
            B = signals.approxAR(Y,p.P);
            c = create_c_vec(p.N,B)*p.gamma;
            A = create_constraint_matrix(p.N);
            b = create_constraint_vec(p.N,p.N);

            options = optimset('Display', 'off');
            [phi,~] = quadprog(Q,c,A(p.N+2:end,:),b(p.N+2:end),A(1:p.N+1,:),b(1:p.N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
            
            l = create_dup_matrix(p.N)*phi;
            L = convert_to_matrix(l);

        end
        function L = AR_mean(Y,p)
            B = zeros(p.N,p.P,p.intervals);
            for j = 1:p.intervals
                [B(:,:,j),~] = signals.approxAR(Y(:,(1+(j-1)*p.interval_length):(j*p.interval_length)),p.P);
            end
            B = mean(B,3);
            B = signals.z_score(B);
            Q = create_Q_matrix(p.N,1);
            c = create_c_vec(p.N,B)*p.gamma;
            A = create_constraint_matrix(p.N);
            b = create_constraint_vec(p.N,p.N);

            options = optimset('Display', 'off');
            [phi,~] = quadprog(Q,c,A(p.N+2:end,:),b(p.N+2:end),A(1:p.N+1,:),b(1:p.N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),zeros(size(Q,1),1),options); %% switched to quadprog...
            l = create_dup_matrix(p.N)*phi;
            L = convert_to_matrix(l);


        end
        function L = AR_mean_L1(Y,p)
            B = zeros(p.N,p.P,p.intervals);
            for j = 1:p.intervals
                [B(:,:,j),~] = signals.approxAR(Y(:,(1+(j-1)*p.interval_length):(j*p.interval_length)),p.P);
            end
            B = mean(B,3);
            B = signals.z_score(B);
            c = create_c_vec_L1(p.N,B,p.gamma);
            Q = zeros(length(c));

            A = create_constraint_matrix(p.N);
            b = create_constraint_vec(p.N,p.N);

            options = optimset('Display', 'off');
            [phi,~] = linprog(c,A(p.N+2:end,:),b(p.N+2:end),A(1:p.N+1,:),b(1:p.N+1),-Inf(size(Q,1),1),Inf(size(Q,1),1),options); %% switched to quadprog...
            l = create_dup_matrix(p.N)*phi;
            L = convert_to_matrix(l);
            % stopped using osqp due to constraints failing


        end
        function [L,Y,L_harvard] = dong(X_noisy,param,beta)
            N = param.N;
            max_iter = param.max_iter;
            alpha = param.alpha;
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
        function [L] = dong_mean(X_noisy,param)
            L_tmp = zeros(param.N,param.N,param.intervals);
            for u = 1:param.intervals
                N = param.N;
                max_iter = param.max_iter;
                alpha = param.alpha;
                beta = param.beta;
                objective = zeros(max_iter,1);
                Y_0 = X_noisy;
                Y = Y_0;
                for i = 1:max_iter
                    L = optimize_laplacian_gaussian(N,Y,alpha,beta);
                    if i == 1
                        L_harvard = L;
                    end
                    R = chol(eye(N) + alpha*L);
                    Y = R \ (R' \ (Y_0));
                    objective(i) = norm(Y-Y_0,'fro')^2 + alpha*vec(Y*Y')'*vec(L) + beta*(norm(L,'fro')^2);
                    if i>=2 && abs(objective(i)-objective(i-1))<10^(-4)
                        break
                    end
                end
                L_tmp(:,:,u) = L;
            end
            L = mean(L_tmp,3);
        end
        function L = threshold(L,t)
            L(abs(L)<t) = 0;
        end
        function t = get_threshold(weights,pred_edges,s)
            sorted_weights = sort(weights,'descend');
            t = sorted_weights(pred_edges);
            if s > .4 && s < .6
                t = t*(1-.5*(s-.4));
            elseif s >= .6
                t = t*.9;
            end
        end
        function g = get_gamma(s)
            a = -0.2724;
            b = 2.6540;
            g = a*exp(b.*s);
            g = 10^g;
        end
        
        function AR_params = create_default_params(signal_params)
            % sets up AR parameters, but gamma and threshold still need to be chosen
            AR_params.N = signal_params.N;
            AR_params.intervals = signal_params.intervals;
            AR_params.interval_length = signal_params.interval_length;
            AR_params.M = signal_params.M;
            AR_params.P = 40;
            AR_params.gamma = 0;
            AR_params.threshold = 0;
        end
        function plot_F(F,N,gammas,thresholds)
            [X_axis,Y_axis] = meshgrid(log10(gammas),thresholds*N);
            s = surf(X_axis,Y_axis,F');
            s.EdgeColor = "none";
            xlabel("log$$_{10}(\gamma)$$",'Interpreter','latex','FontSize',14)
            ylabel("$$\tilde{\delta}$$",'Interpreter','latex','FontSize',14)
            zlabel("F-score",'Interpreter','latex','FontSize',14)
        end
        function [fits,P] = sweep_fits(y,Pmax)
            P = 4:4:Pmax;
            fits = zeros(length(P),1);
            for i = 1:length(P)
                [~,fits(i)] = signals.approxAR(y,P(i));
            end
        end
    end
end

