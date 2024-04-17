function F = calculate_F_score(theoretical, experimental)
    P = calculate_precision(theoretical,experimental);
    R = calculate_recall(theoretical,experimental);
    F = 2/(P^-1+R^-1);
end