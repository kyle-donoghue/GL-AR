function R = calculate_recall(theoretical,experimental)
    [tp,~,~,fn] = calculate_result_true_false(theoretical,experimental);
    R = tp/(tp+fn);
end