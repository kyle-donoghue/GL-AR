function signalPlot(x)
    N = size(x,2);
    h = plot(x);
    leg=repmat({'Y'},N,1);
    leg=genvarname(leg);
    legend(leg);
    for i = 1:numel(h)
        % Add a new row to DataTip showing the DisplayName of the line
        h(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('',repmat({h(i).DisplayName},size(h(i).XData))); 
    end
end