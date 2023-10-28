close all hidden;
[X_axis,Y_axis] = meshgrid(log(gammas),thresholds);
F2 = permute(F,[1 3 2]);
figure;
for i = [1:length(noises)]
    s = surf(X_axis,Y_axis,F2(:,:,i)');
    s.EdgeColor = 'none';
    [noises(i) max(F2(:,:,i),[],"all")]
    pause(2)
end
%%
[maxF,ind] = max(F2,[],"all")
[row,col,page] = ind2sub(size(F2),ind)
