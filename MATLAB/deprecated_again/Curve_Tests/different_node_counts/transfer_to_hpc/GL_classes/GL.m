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
        function L = AR_mean(Y,p)
            B = zeros(p.N,p.P,p.intervals);
            for j = 1:p.intervals
                B(:,:,j) = signals.approxAR(Y(:,(1+(j-1)*p.interval_length):(j*p.interval_length)),p.P);
            end
            B = mean(B,3);
            B = signals.z_score(B);
            Q = create_Q_matrix(p.N,1);
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
                L_tmp(:,:,u) = L;
            end
            L = mean(L_tmp,3);
        end
        function L = threshold(L,t)
            L(abs(L)<=t) = 0;
        end
        function t = get_threshold(edge_spread,N)
            % a = 21.34;
            % b = -13.88; % switching because it seemed off 12/19

            % a =       17.78;
            % b =      -14.13;
            % t_offset = .75;
            % t = (a*exp(b*(edge_spread/N))+t_offset)/N; % switching to exp2

            % a =       4.225;
            % b =      -4.762;
            % c =       19.85;
            % d =      -32.88;
            % t_offset = .25;
            % 
            % c = c*20/N;
            % d = d*sqrt(N/20);
            % 
            % t = (a*exp(b*(edge_spread/N))+c*exp(d*(edge_spread/N))+t_offset)/N;
            
            %switching to exp1 with cliff constraints using edge_spread2
            % a =       5.241;
            % b =      -1.264;
            % a = 5.979;
            % b = -1.225;
            % t_offset = .75;
            % t = (a*exp(b*(edge_spread/N))+t_offset)/N;

            %switching to exp2 again but using cliff constraints and edge_spread2
            % a = 16.1429;
            % b = -5.9860;
            % c = 4.1450;
            % d = -1.1352;
            % t_offset = .8;
            % t = (a*exp(b*(edge_spread/N))+c*exp(d*(edge_spread/N))+t_offset)/N;

            % de-normalizing edge_spread2 for better N scaling
            % adding N scaling for different threshold curves
            % a = 16.1429*(N/20);
            % b = -0.2993*(20/N)^1;
            % c = 4.1450*(N/20);
            % d = -0.0568*(20/N)^1;
            % % t_offset = .8;
            % t_offset = .9;
            % t = (a*exp(b*(edge_spread))+c*exp(d*(edge_spread))+t_offset)/N;

            % a = 22.24*(N/20);
            % b = -0.5748*(20/N)^1.75;
            % c = 6.243*(N/20);
            % d =  -0.07697*(20/N)^1;
            % t_offset = .95;
            % t = (a*exp(b*(edge_spread))+c*exp(d*(edge_spread))+t_offset)/N;
            
            a = 1.1120*(N);
            b = -108.7225/(N^1.75);
            c = 0.3122*(N);
            d =  -1.5394/N;
            t_offset = .95;
            t = (a*exp(b*(edge_spread))+c*exp(d*(edge_spread))+t_offset)/N;

        end
        function g = get_gamma(s,t,N,v_std)
            cliff_curve_offset = 1;
            cliff_curve_delay = min(2.25,max(.775,v_std));
            ridge_offset = .9;
            t_offset = .95;

            [a,b] = GL.get_cliff_curve(s);
            % b=b*max(1,N/20);
            
            t = t*N; % get normalized threshold, since cliff curve is normalized
            % cutoff_multiplier = 1+(.0005*(N/20));
            cutoff_multiplier = 1+(2.5e-05*(N));
            if t > t_offset*cutoff_multiplier
                g = max((log((t-ridge_offset*cliff_curve_offset)/a)/b+cliff_curve_delay),0);
            else
                g = 0;
            end
            g = 10^g;
        end
        function [a,b] = get_cliff_curve(s)
            % a_a = 0.1289;
            % a_b = -5.855;
            % a_offset = 0.0060;
            % b_a1 =  0.4011;
            % b_c1 = 0.2049;
            % b_offset = 1.5;

            % below is F_t=.02
            % a_a = 2.579;
            % a_b = -5.855;
            % a_offset = 0.1190;
            % b_a1 =  0.3925;
            % b_c1 = 0.2044;
            % b_offset = 1.5;

            % below is F_t=.001
            a_a = 2.665;
            a_b = -4.314;
            a_offset = 0.1190;
            b_a1 =  0.5195;
            b_c1 = 0.2875;
            b_offset = 1.45;

            a = a_a*exp(a_b*s)+a_offset;
            b = b_a1*exp(-(s/b_c1)^2)+b_offset;
        end
        function [t1,t2] = plot_cliff_curve(g,a,b,N,v_std)
            % b=b*max(1,N/20);
            ridge_offset = .9;
            cliff_curve_offset = 1;
            cliff_curve_delay = min(2,max(.775,v_std));
            t1 = a*exp(b*(log10(g)-cliff_curve_delay))+ridge_offset*cliff_curve_offset;
            t2 = a*exp(b*log10(g))+cliff_curve_offset;
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
    end
end

