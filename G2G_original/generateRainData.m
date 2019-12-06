function rainfallRate = generateRainData(avgSunnyPeriod, mu, sigma, nTimeSteps)
    rainfallRate = zeros(nTimeSteps, 1);
    isRain = false;
    for time = 1:nTimeSteps
        if isRain
            rainfallRate(time) = rainfallRate(time - 1) + normrnd(-mu, sigma);
            if rainfallRate(time) < 0
                rainfallRate(time) = 0;
                isRain = false;
            end
        elseif rand < 1/avgSunnyPeriod
                isRain = true;
        end
    end
end

