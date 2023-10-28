close all hidden;

for k = 1:length(gammas)
    [max_val,ind] = max(Fmeasure(:,:,k),[],[1 2],"linear");
    [r,c] = ind2sub(size(Fmeasure(:,:,k)),ind);
    f_gamma(k) = max_val;
    p_gamma(k) = precision(r,c,k);
    r_gamma(k) = recall(r,c,k);
end
figure;hold on;
plot(gammas,f_gamma);
plot(gammas,p_gamma);
plot(gammas,r_gamma);

for k = 1:length(gammas)
    p_gamma(k) = max(precision(:,:,k),[],[1 2]);
    r_gamma(k) = max(recall(:,:,k),[],[1 2]);
    f_gamma(k) = 2/(p_gamma(k)^-1+r_gamma(k)^-1);
end
figure;hold on;
plot(gammas,f_gamma);
plot(gammas,p_gamma);
plot(gammas,r_gamma);


figure;
for i = [1:length(gammas)]
    surf(X_axis,Y_axis,Fmeasure(:,:,i)')
    [gammas(i) max(Fmeasure(:,:,i),[],"all")]
    pause(.5)
end
%%
[maxF,ind] = max(Fmeasure,[],"all")
[row,col,page] = ind2sub(size(Fmeasure),ind)
maxGamma = gammas(page)
maxThresh = thresholds(row)
maxAvg = thresholds_avg(col)