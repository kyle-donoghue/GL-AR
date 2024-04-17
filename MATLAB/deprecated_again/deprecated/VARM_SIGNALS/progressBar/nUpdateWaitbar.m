function nUpdateWaitbar(N)
    persistent p g loading
    if isempty(p)
        p = 1;
    end
    if isempty(g) || ~isvalid(g)
        g = waitbar(0, 'Please wait ...');
        loading = tic;
    end
    ticktime = toc(loading)/60;
    left = ((N - p)*ticktime);
    waitbar(p/(N), g, sprintf('%12.9f , %12f',ticktime,left));
    loading = tic;
    p = p + 1;
end