classdef signals
    %SIGNAL Various methods to create and approximate signals

    methods (Static=true)
        function B = approxAR(x,P)
            N = size(x,1);
            B = zeros(N,P+1);
            for i = 1:N
                B(i,:) = ar(x(i,:),P,'fb').A; % switched to this for faster speed
            end
            B = B(:,2:end);
        end
        function [V,I,act,freqs] = circuit_sine(p,L)
            act = sort(randperm(p.N,p.active));
            freqs = generate_frequencies2(p.active,p.min_freq,p.max_freq);
            if p.isComplex
                compl_I = zeros(p.N,p.active);
                for i = 1:p.active
                    compl_I(act(i),i) = 1;
                end
                compl_v = pinv(L)*compl_I;
                V = zeros(p.N,p.M);
                I = zeros(p.N,p.M);
                I(act,:) = sin(2*pi*freqs*(1/p.Fs:1/p.Fs:p.M/p.Fs));
                for i = 1:p.N %doesn't work for activated nodes > 1
                    for j = 1:p.active
                        V(i,:) = V(i,:) + abs(compl_v(i,j))*sin(2*pi*freqs(j)*(1/p.Fs:1/p.Fs:p.M/p.Fs)+angle(compl_v(i,j)));
                    end
                end
            else
                I = zeros(p.N,p.M);
                I(act,:) = sin(2*pi*freqs*(1/p.Fs:1/p.Fs:p.M/p.Fs));
                V = pinv(L)*I;
            end
            V = signals.add_noise(p,V);
        end
        function x = z_score(x) %expects NxM signal to z_score across rows
            x = (x-mean(x,2))./std(x,0,2);
        end
        function plot(x)
            N = size(x,2);
            h = plot(x);
            leg = strcat(repmat('Y',N,1),num2str((1:N)'));
            legend(leg);
            for i = 1:numel(h)
                % Add a new row to DataTip showing the DisplayName of the line
                h(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('',repmat({h(i).DisplayName},size(h(i).XData))); 
            end
        end
        function Y = filterTensorSpectrum(p,X,G)
            % tensor spectrum X, and graph tensor G, to make tensor spectrum Y
            Y = zeros(p.N,p.intervals,p.interval_length);
            for i = 1:p.interval_length
                Y(:,:,i) = G(:,:,i)*X(:,:,i);
            end
        end
        function X = createTensorSpectrum(p,x)
            % create tensor spectrum from x
            X = zeros(p.N,p.intervals,p.interval_length);
            for i = 1:p.intervals
                tmp = fft(x(:,(1+(i-1)*p.interval_length):(i*p.interval_length)).');
                X(:,i,:) = tmp.';
            end
        end
        function x = inverseTensorSpectrum(p,X)
            x = zeros(p.N,p.interval_length,p.intervals);
            for i = 1:p.intervals
                tmp = permute(X(:,i,:),[3 1 2]);
                x(:,:,i) = real(ifft(tmp)).';
            end
            x = reshape(x,p.N,p.M);
        end
        function x = add_noise(p,x)
            x = x+var(x,0,2)/p.SNR.*randn(p.N,p.M);
        end
        function [fits,aics,p] = test_fit(x)
            p = 2:2:floor(length(x)/2);
            [fits,aics] = deal(zeros(length(p),1));
            for i = 1:length(p)
                report = ar(x,p(i),'fb').Report.Fit;
                fits(i) = report.FitPercent;
                aics(i) = report.AIC;
            end
        end
        function x = create_raw_sine(p)
            freqs = generate_frequencies2(p.N,p.min_freq,p.max_freq);
            x = sin(2*pi*freqs*(1/p.Fs:1/p.Fs:p.M/p.Fs));
        end
        function a = create_IIR(p)
            pole_list = max(.45, min(.95, (p.pole_variance*randn(p.order,1)+p.pole_mean))).*exp(1j*pi*rand(p.order,1)); % saturate pole to +/- .95
            pole_list2 = [pole_list;pole_list];
            pole_list2(length(pole_list)+1:end) = conj(pole_list2(length(pole_list)+1:end));
            a = poly(pole_list2);
            H = freqz(1,a);
            a = a*max(abs(H));
            % b = 1/max(abs(H));
            % a = [b a];
        end
        function y = customFilter(p,x,D)
            order = size(D,3);
            y = zeros(p.N,p.M);
            for i = 1:p.M
                y(:,i) = x(:,i);
                for k = 1:order-1
                    index = i-(p.delay+(k-1));
                    if index >= 1
                        if k == 1
                            y(:,i) = y(:,i) + D(:,:,k)*x(:,index);
                        else
                            y(:,i) = y(:,i) + D(:,:,k)*y(:,index);
                        end
                    end
                end
            end
        end
        function x = createFullSinePulse(p)
            N = p.N;
            M = p.M;
            Fs = p.Fs;
            freqs = generate_frequencies2(N,p.min_freq,p.max_freq);
            minSeparation = p.minSeparation;
            maxSeparation = p.maxSeparation;

            x = zeros(N,M);
            for i = 1:N
                inSignal = 1;
                f = freqs(i);
                indices = [randi(maxSeparation)];
                while inSignal
                   indices(end+1) = indices(end)+floor(Fs/f)+minSeparation+randi(maxSeparation);
                   inSignal = (indices(end)+floor(Fs/f)) < M;
                end
                indices = indices(1:end-1);
                for k = 1:length(indices)
                    x_tmp = signals.createPulse(f,Fs);
                    x(i,indices(k):indices(k)+length(x_tmp)-1) = x_tmp;
                end
            end
        end
        function x = createPulse(freq,Fs)
            N = floor(Fs/freq);
            t = 0:1/Fs:(N-1)/Fs;
            x = sin(2*pi*freq*t);
        end
    end
end

