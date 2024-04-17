close all hidden;


figure;
for i = [1:length(gammas)]
    s = surf(X_axis,Y_axis,F(:,:,i)');
    s.EdgeColor = 'none';
    [shunts(i) max(F(:,:,i),[],"all")]
    pause(.5)
end
%%
[maxF,ind] = max(Fmeasure,[],"all")
[row,col,page] = ind2sub(size(Fmeasure),ind)
maxGamma = gammas(page)
maxThresh = thresholds(row)
maxAvg = thresholds_avg(col)