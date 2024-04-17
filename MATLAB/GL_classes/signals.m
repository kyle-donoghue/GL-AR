classdef signals
    %SIGNAL Various methods to create and approximate signals

    methods (Static=true)
        function [B,fit] = approxAR(x,P)
            N = size(x,1);
            B = zeros(N,P+1);
            fit = zeros(N,1);
            for i = 1:N
                mdl = ar(x(i,:),P,'fb');
                B(i,:) = mdl.A; % switched to this for faster speed
                fit(i) = mdl.Report.Fit.nAIC;
            end
            B = B(:,2:end);
            fit = mean(fit);
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
                for i = 1:p.N 
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
        function [x,mu,sigma] = z_score(x) %expects NxM signal to z_score across rows
            mu = mean(x,2);
            sigma = std(x,0,2);
            x = (x-mu)./sigma;
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
        function X_noisy = createGFTsignals(p,L_0)
            [V,D] = eig(full(L_0));
            sigma = pinv(D);
            mu = zeros(1,p.N);
            gftcoeff = mvnrnd(mu,sigma,p.M);
            X = V*gftcoeff';
            X_noisy = signals.add_noise(p,X);
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
            variance = mean(var(x,0,2));
            x = x+variance/p.SNR.*randn(p.N,p.M);
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
            pole_list = max(.2, min(.8, (p.pole_variance*randn(p.order,1)+p.pole_mean))).*exp(1j*pi*rand(p.order,1)); % saturate pole to +/- .95
            pole_list2 = [pole_list;pole_list];
            pole_list2(length(pole_list)+1:end) = conj(pole_list2(length(pole_list)+1:end));
            a = poly(pole_list2);
            H = freqz(1,a);
            a = a*max(abs(H));
        end
        function x = createARprocesses(p)
            p.order = p.processOrder;
            x = randn(p.N,p.M);
            for i = 1:p.N
                a = signals.create_IIR(p);
                x(i,:) = filter(1,a,x(i,:));
            end
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
        function x = createFullSinePulse(p,pulses)
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
                   indices(end+1) = indices(end)+pulses*floor(Fs/f)+minSeparation+randi(maxSeparation);
                   inSignal = (indices(end)+floor(Fs/f)) < M;
                end
                indices = indices(1:end-1);
                for k = 1:length(indices)
                    x_tmp = signals.createPulse(f,Fs);
                    x(i,indices(k):indices(k)+pulses*length(x_tmp)-1) = repmat(x_tmp,1,pulses);
                end
            end
        end
        function x = createFullNyquistPulse(p)
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
                   indices(end+1) = indices(end)+floor(4*Fs/f)+minSeparation+randi(maxSeparation);
                   inSignal = (indices(end)+floor(4*Fs/f)) < M;
                end
                indices = indices(1:end-1);
                for k = 1:length(indices)
                    x_tmp = signals.createNyquistPulse(f,Fs);
                    x(i,indices(k):indices(k)+length(x_tmp)-1) = x_tmp;
                end
            end
        end
        function x = createFullSinePulseWindow(p)
            N = p.N;
            M = p.M;
            Fs = p.Fs;
            freqs = generate_frequencies2(N,p.min_freq,p.max_freq);
            minSeparation = p.minSeparation;
            maxSeparation = p.maxSeparation;
            minPer = p.minPer;
            maxPer = p.maxPer;
            
            x = zeros(N,M);
            for i = 1:N
                inSignal = 1;
                f = freqs(i);
                indices = [randi(maxSeparation)];
                currPer = [];
                while inSignal
                   currPer(end+1) = minPer + (maxPer-minPer).*rand(1);
                   indices(end+1) = indices(end)+floor(Fs*currPer(end)/f)+minSeparation+randi(maxSeparation);
                   inSignal = (indices(end)+floor(Fs*currPer(end)/f)) < M;
                end
                indices = indices(1:end-1);
                for k = 1:length(indices)
                    x_tmp = signals.createWindowPulse(f,Fs,currPer(k));
                    x(i,indices(k):indices(k)+length(x_tmp)-1) = x_tmp;
                end
            end
        end
        function x = createFullRectPulse(p)
            N = p.N;
            M = p.M;
            Fs = p.Fs;
            % periods = 10+(1:p.N);
            periods = mod(10+(1:p.N),30);
            minSeparation = p.minSeparation;
            maxSeparation = p.maxSeparation;

            x = zeros(N,M);
            for i = 1:N
                inSignal = 1;
                period = periods(i);
                indices = [randi(maxSeparation)];
                while inSignal
                   indices(end+1) = indices(end)+period+minSeparation+randi(maxSeparation);
                   inSignal = (indices(end)+period) < M;
                end
                indices = indices(1:end-1);
                for k = 1:length(indices)
                    x_tmp = ones(1,period);
                    x(i,indices(k):indices(k)+length(x_tmp)-1) = x_tmp;
                end
            end
        end
        function x = createFullRectPulseWindow(p)
            N = p.N;
            M = p.M;
            Fs = p.Fs;
            periods = 50+(1:p.N);
            minSeparation = p.minSeparation;
            maxSeparation = p.maxSeparation;

            x = zeros(N,M);
            for i = 1:N
                inSignal = 1;
                period = periods(i);
                indices = [randi(maxSeparation)];
                while inSignal
                   indices(end+1) = indices(end)+period+minSeparation+randi(maxSeparation);
                   inSignal = (indices(end)+period) < M;
                end
                indices = indices(1:end-1);
                for k = 1:length(indices)
                    x_tmp = ones(1,period).*tukeywin(period)';
                    x(i,indices(k):indices(k)+length(x_tmp)-1) = x_tmp;
                end
            end
        end
        function y_noisy = generateFilteredRectPulse(p,G)
            x = signals.createFullRectPulse(p);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredRectPulseWindow(p,G)
            x = signals.createFullRectPulseWindow(p);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredSinePulse(p,G,pulses)
            x = signals.createFullSinePulse(p,pulses);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredSinePulseWindow(p,G)
            x = signals.createFullSinePulseWindow(p);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredSine(p,G)
            x = signals.create_raw_sine(p);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredNyquistPulse(p,G)
            x = signals.createFullNyquistPulse(p);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredARProcess(p,G)
            x = signals.createARprocesses(p);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function y_noisy = generateFilteredGaussian(p,G)
            x = randn(p.N,p.M);
            X = signals.createTensorSpectrum(p,x);
            Y = signals.filterTensorSpectrum(p,X,G);
            y = signals.inverseTensorSpectrum(p,Y);
            y_noisy = signals.add_noise(p,y);
        end
        function x = createPulse(freq,Fs)
            N = floor(Fs/freq);
            t = 0:1/Fs:(N-1)/Fs;
            x = sin(2*pi*freq*t);
        end
        function x = createNyquistPulse(freq,Fs)
            N = floor(4*Fs/freq);
            t = -(N/Fs/2):1/Fs:(N-1)/Fs/2;
            x = sinc(2*freq*(t));
            x = x.*gausswin(N)';
        end
        function x = createWindowPulse(freq,Fs,per)
            N = floor(Fs*per/freq);
            t = 0:1/Fs:(N-1)/Fs;
            x = sin(2*pi*freq*t);
            x = x.*tukeywin(N)';
        end
        function signal_params = create_default(N,graph_type)
            signal_params.N = N;
            signal_params.raw = 1;
            signal_params.directed = 0;
            
            signal_params.SNR = 2e0;
            signal_params.Fs = 100;
            signal_params.min_freq = 2;
            signal_params.max_freq = signal_params.Fs/2;
            
            signal_params.intervals = min(20,signal_params.N);
            % signal_params.intervals = signal_params.N;
            signal_params.interval_length = 512;
            signal_params.M = signal_params.intervals*signal_params.interval_length;
            
            signal_params.order = 1;
            signal_params.zero_variance = 1;
            signal_params.zero_mean = 1;
            signal_params.pole_variance = .1;
            signal_params.pole_mean = .5;

            signal_params.processOrder = 3;
            
            signal_params.min_delay = 5;
            signal_params.max_delay = 15;
            
            
            signal_params.minSeparation = 50;
            signal_params.maxSeparation = 200;

            signal_params.minPer = 20;
            signal_params.maxPer = 50;

            signal_params.graph_type = graph_type;
        end
        function p = create_empty(sz)
            p.N = sz(1);
            p.M = sz(2);
            p.intervals = 1;
            p.interval_length = p.M;
        end
    end
end

