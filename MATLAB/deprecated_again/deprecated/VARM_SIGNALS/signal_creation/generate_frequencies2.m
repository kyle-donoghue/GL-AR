function p = generate_frequencies2(numFrequencies,minFrequency,maxFrequency)
    p = primes(maxFrequency*10);
    p = p(p>minFrequency*10);
    p = p(randperm(length(p),numFrequencies))'/10;
    p = sort(p);
end