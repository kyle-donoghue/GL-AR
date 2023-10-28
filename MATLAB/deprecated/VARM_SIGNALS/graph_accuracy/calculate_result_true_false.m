function [tp,tn,fp,fn] = calculate_result_true_false(theoretical, experimental)
    tp = 0;
    tn = 0;
    fp = 0;
    fn = 0;
    n = size(theoretical,1);
    for i = 1:n
        for j = 1:i-1 % iterate through lower triangle
            error = 2*double(theoretical(i,j)>0) - double(experimental(i,j)>0);
            % tp = 2-1 = 1,
            % tn = 0-0 = 0,
            % fp = 0-1 = -1,
            % fn = 2-0 = 2
            %case statement  for accuracy
            switch error
                case 1
                    tp = tp+1;
                case 0
                    tn = tn+1;
                case -1
                    fp = fp+1;
                case 2
                    fn = fn+1;
            end
        end
    end
end