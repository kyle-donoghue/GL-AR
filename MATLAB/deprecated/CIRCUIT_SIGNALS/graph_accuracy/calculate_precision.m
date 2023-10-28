function P = calculate_precision(theoretical,experimental)
    [tp,~,fp,~] = calculate_result_true_false(theoretical,experimental);
    P = tp/(tp+fp);
end