function nonHarmonicFrequencies = generate_frequencies(numFrequencies,minFrequency,maxFrequency,tolerance)
    % Initialize an empty array for non-harmonic frequencies
    nonHarmonicFrequencies = [];
    % Loop until you have the desired number of non-harmonic frequencies
    while length(nonHarmonicFrequencies) < numFrequencies
        % Generate a random frequency within the specified range
        randomFrequency = minFrequency + (maxFrequency - minFrequency) * rand();
        
        % Check if the random frequency is a harmonic of any existing frequency
        isHarmonic = false;
        for j = 1:length(nonHarmonicFrequencies)
            % Check if the current frequency is within a tolerance of being a harmonic
            if abs(randomFrequency / nonHarmonicFrequencies(j) - round(randomFrequency / nonHarmonicFrequencies(j))) < tolerance
                isHarmonic = true;
                break;
            end
        end
        
        % If it's not a harmonic, add it to the list of non-harmonic frequencies
        if ~isHarmonic
            nonHarmonicFrequencies = [nonHarmonicFrequencies;randomFrequency];
        end
    end
    nonHarmonicFrequencies = sort(nonHarmonicFrequencies);
end