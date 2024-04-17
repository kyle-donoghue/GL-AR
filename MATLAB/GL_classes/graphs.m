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
            heatmap(abs(L),'CellLabelColor','none');
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
            [nonzeros_row, nonzeros_col] = find(triu(A)>0);
            for i = 1:length(nonzeros_row)
                a = signals.create_IIR(p);
                b = [zeros(1,randi([p.min_delay p.max_delay],1)) 1];
                [H,w] = freqz(b,a,p.interval_length,'whole');
                
               
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
        function s = get_density(graph_type,graph_param,N)
            if strcmp(graph_type,'er')
                s = graph_param;
            elseif strcmp(graph_type,'pa')
                s = 2*graph_param/N;
            elseif strcmp(graph_type,'gaussian')
                p1 = -0.7808;
                p2 = 1.535;
                p3 = 0.05708;
                s = p1*graph_param.^3+p2*graph_param.^2+p3*graph_param;
            end
        end
        function n = unconnected_vertices(L)
            n = sum(diag(L) == 0);
        end
        function weights = get_weights(L)
            N = length(L);
            weights = zeros(N*(N-1)/2,1);
            
            count = 1;
            for row = 1:(N-1)
                for col = row+1:N
                    weights(count) = L(row,col);
                    count = count+1;
                end
            end
            weights = -weights;
        end
        function m = max_edges(N)
            m = (N*(N-1)/2);
        end
        function d = density(L)
            w = graphs.get_weights(L);
            d = nnz(w)/length(w);
        end
        function F = minFscore(density)
            F = 2*(density)/(1+density);
        end
    end
end

