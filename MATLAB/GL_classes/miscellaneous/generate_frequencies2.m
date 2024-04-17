function p = generate_frequencies2(numFrequencies,minFrequency,maxFrequency)
    scalar = 25;
    p = primes(maxFrequency*scalar);
    p = p(p>minFrequency*scalar);
    p = p(randperm(length(p),numFrequencies))'/scalar;
    p = sort(p);
end