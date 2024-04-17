classdef graphs
    %GRAPH methods to create and compare graphs
    
    methods (Static=true)
        function [L,A,A_d] = create(p,varargin1)
            opt = p.graph_type;
            if ~isfield(p,'raw')
                raw = 0;
            else
                raw = p.raw;
            end
            if ~isfield(p,'directed')
                directed = 0;
            else
                directed = p.directed;
            end
            N = p.N;
            if ~raw
                trace_norm = p.trace_normalization;
                diag_add = p.additional_diag;
                isComplex = p.isComplex;
            end
            % if nargin == 1
            %     if strcmp(opt,'chain') == 0
            %         error('number of input variables not correct :(')
            %     end
            % elseif nargin == 2
            %     if strcmp(opt,'gaussian') || strcmp(opt,'ff')
            %         error('number of input variables not correct :(')
            %     end
            % elseif nargin == 3
            %     if strcmp(opt,'er') || strcmp(opt,'pa')
            %         error('number of input variables not correct :(')
            %     end
            % end
            %% generate coordinates of vertices
            plane_dim = 1;
            XCoords = plane_dim*rand(N,1);
            YCoords = plane_dim*rand(N,1);
            %% construct the graph
            switch opt
                case 'gaussian' % random graph with Gaussian weights
                    % T = varargin1;
                    T = 0.75;
                    % s = varargin2; %making it easier for default graphs
                    s = varargin1;
                    d = distanz([XCoords,YCoords]'); 
                    W = exp(-d.^2/(2*s^2)); 
                    W(W<T) = 0; % Thresholding to have sparse matrix
                    W = 0.5*(W+W');
                    G = W-diag(diag(W));
                    
                case 'er' % Erdos-Renyi random graph
                    p_er = varargin1;
                    G = erdos_reyni(N,p_er);
                    
                case 'pa' % scale-free graph with preferential attachment
                    m = varargin1;
                    G = preferential_attachment_graph(N,m);
                    
                % case 'ff', % forest-fire model
                %     p = varargin1;
                %     r = varargin2;
                %     G = forest_fire_graph(N,p,r);
                    
                case 'chain' % chain graph
                    G = spdiags(ones(N-1,1),-1,N,N);
                    G = G + G';
            end
            A = full(G);
            if ~raw
                if isComplex
                    phasors = 2*pi.*rand(N);
                    A_tmp = triu(A).*exp(1j*phasors);
                    A = A_tmp.'+triu(A_tmp,1);
                end
            end
            L = graphs.to_laplacian(A);
            % L = diag(sum(abs(A)))-A;
            if ~raw
                L = L*trace_norm/trace(L);
                L = L+diag_add*eye(N);
            end
            if directed
                A_d = graphs.to_directed(A);
            else
                A_d = A;
            end
            while (~any(L,'all') || all(L,'all'))
                [L,A,A_d] = graphs.create(p,varargin1);
            end

        end
        function [precision,recall,f,NMI,num_of_edges] = performance(L_0,L)            
            L_0tmp = L_0-diag(diag(L_0));
            edges_groundtruth = squareform(L_0tmp)~=0;
            Ltmp = L-diag(diag(L));
            edges_learned = squareform(Ltmp)~=0;
            num_of_edges = sum(edges_learned);
            if num_of_edges > 0
                [precision,recall] = perfcurve(double(edges_groundtruth),double(edges_learned),1,'Tvals',1,'xCrit','prec','yCrit','reca');
                if precision == 0 && recall == 0
                    f = 0;
                else
                    f = 2*precision*recall/(precision+recall);
                end
                NMI = perfeval_clus_nmi(double(edges_groundtruth),double(edges_learned));
            else
                precision = 0;
                recall = 0;
                f = 0;
                NMI = 0;
            end
        end
        function num_of_edges = num_edges(L)
            Ltmp = L-diag(diag(L));
            edges_learned = squareform(Ltmp)~=0;
            num_of_edges = sum(edges_learned);
        end
        function plot(L,t)
            heatmap(abs(L));
            colormap parula;
            title(t);
        end
        function vcompare(L1,L2,t1,t2)
            subplot(1,2,1);
            graphs.plot(L1,t1);
            subplot(1,2,2);
            graphs.plot(L2,t2);
        end
        function A = to_adjacency(L)
            A = diag(diag(L))-L;
        end
        function L = to_laplacian(A)
            L = diag(sum(A)) - A;
        end
        function A = to_directed(A)
            N = size(A,2);
            for i = 2:N
                for j = 1:i-1
                    if rand >= .5
                        A(i,j) = 0;
                    else
                        A(j,i) = 0;
                    end
                end
            end
        end
        function G = createGraphTensor(p,A)
            G = zeros([size(A) p.interval_length]);
            [nonzeros_row, nonzeros_col] = find(A>0);
            for i = 1:length(nonzeros_row)
                % zero_list = (p.zero_variance*randn(p.order,1)+p.zero_mean).*exp(0*1j*2*pi*rand(p.order,1));
                % pole_list = max(-.95, min(.95, (p.pole_variance*randn(p.order,1)+p.pole_mean))).*exp(0*1j*2*pi*rand(p.order,1)); % saturate pole to +/- .95
                % b = poly(zero_list);
                % a = poly(pole_list);
                a = signals.create_IIR(p);
                % b = 1;
                % b = [zeros(1,30) 1];
                b = [zeros(1,randi([p.min_delay p.max_delay],1)) 1];
                H = freqz(b,a,p.interval_length,'whole');
                
                % normalize to 0->1
                % H = H/(max(H)-min(H));
                % H = H-min(H);
               
                G(nonzeros_row(i),nonzeros_col(i),:) = H;
                G(nonzeros_col(i),nonzeros_row(i),:) = H; % this line is for undirected
       
            end
            for i = 1:p.N % preserve signal on itself
                G(i,i,:) = ones(p.interval_length,1);
            end
        end
        function D = create_VARM_tensor(p,A)
            D = zeros([size(A) p.order*2+2]);
            [nonzeros_row, nonzeros_col] = find(A>0);
            for i = 1:length(nonzeros_row)
                a = signals.create_IIR(p);
                D(nonzeros_row(i),nonzeros_col(i),:) = a;
            end
        end
        function view_tensor_spectrum(G)
            figure;hold on;
            [nonzeros_row, nonzeros_col] = find(G(:,:,1)>0);
            for i = 1:length(nonzeros_row)
                if nonzeros_col(i) ~= nonzeros_row(i)
                    H = permute(G(nonzeros_row(i),nonzeros_col(i),:),[3 1 2]);
                    plot(abs(H));
                end
            end
        end
        function s = get_sparsity(graph_type,graph_param)
            if strcmp(graph_type,'er')
                s = graph_param;
            elseif strcmp(graph_type,'pa')
                p1 = -0.002632;
                p2 = 0.1026;
                s = p1*graph_param.^2+p2*graph_param;
            elseif strcmp(graph_type,'gaussian')
                p1 = -0.7808;
                p2 = 1.535;
                p3 = 0.05708;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param;
            end
        end
        function s = get_edge_spread(graph_type,graph_param,N)
            % this is avg/max+avg
            if strcmp(graph_type,'er')
                s = graph_param*N;
            elseif strcmp(graph_type,'pa')
                p1 = -0.002655;
                p2 = 0.1021;
                p3 = 0.00955;
                s = p1*graph_param.^2+p2*graph_param+p3;
                s = s*N;
            elseif strcmp(graph_type,'gaussian')
                p1 = -0.5531;
                p2 = 1.117;
                p3 = 0.1664;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param;
                s = s*N;
            end
        end
        function s = get_edge_spread2(graph_type,graph_param,N)
            % this is avg*max
            % this is based on N=20
            % scaled by N^2 for ER
            % scaled by N^2.5 for RND
            % scaled by N for BA
            if strcmp(graph_type,'er')
                % p1 = -245;
                % p2 = 597.4;
                % p3 = 15.39;
                % p4 = 2.035;
                % s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                % s = s*(N^2/20^2);
                p1 = -0.6125;
                p2 = 1.4935;
                p3 = 0.0385;
                p4 = 0.0051;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                s = s*(N^2);
            elseif strcmp(graph_type,'pa')
                % p1 = -0.1316;
                % p2 = 2.213;
                % p3 = 18.36;
                % p4 = -5.745;
                % s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                % s = s*(N/20);
                p1 = -0.0066;
                p2 = 0.1106;
                p3 = 0.9180;
                p4 = -0.2873;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                s = s*(N);
            elseif strcmp(graph_type,'gaussian')
                % a1 = 278.5;
                % b1 = 1.184;
                % c1 = 0.5191;
                % s = a1*exp(-((graph_param-b1)/c1).^2)-2;
                % s = s*(N^2.5/20^2.5);
                a1 = 0.1557;
                b1 = 1.184;
                c1 = 0.5191;
                s = a1*exp(-((graph_param-b1)/c1).^2)-0.0011;
                s = s*(N^2.5);
            end
        end
        function n = unconnected_vertices(L)
            n = sum(diag(L) == 0);
        end
        function s = get_vertice_std(graph_type,graph_param,N)
            if strcmp(graph_type,'er')
                p1 = 0.437;
                p2 = -6.362;
                p3 = 5.971;
                p4 = 0.7074;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                % s = s*sqrt(20/N);
            elseif strcmp(graph_type,'pa')
                p1 = 0.003242;
                p2 = -0.1035;
                p3 = 0.8698;
                p4 = 1.035;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                % s = s*sqrt(20/N);
            elseif strcmp(graph_type,'gaussian')
                % p1 = -4.656;
                % p2 = 3.654;
                % p3 = 3.577;
                % p4 = 0.1004;
                % s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                % s = s*(N/20);
                p1 = -0.2328;
                p2 = 0.1827;
                p3 = 0.1789;
                p4 = 0.0050;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param+p4;
                s = s*(N);
            end
        end
    end
end

