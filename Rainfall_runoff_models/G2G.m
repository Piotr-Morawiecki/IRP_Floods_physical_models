function [outputArg1,outputArg2] = G2G(domainSize, lengthStep, timeStep, ...
    nTimeSteps, alpha, fastFlowTime, slowFlowTime, maxStorage, rainfallRate)
    % types of stores
    mainStore = zeros(domainSize,1);
    fastStore = zeros(domainSize,1);
    slowStore = zeros(domainSize,1);

    % memory
    mainStoreData = zeros(nTimeSteps, domainSize);
    fastStoreData = zeros(nTimeSteps, domainSize);
    slowStoreData = zeros(nTimeSteps, domainSize);
    dischargeData = zeros(nTimeSteps, 1);

    for time = 1:nTimeSteps
        mainStore = mainStore + timeStep * rainfallRate;
        runoff = max(0, mainStore - maxStorage);
        mainStore = min(maxStorage, mainStore);
        flow = [0; fastStore - [fastStore(2:end); 0]];
        fastStore = fastStore - timeStep * diff(flow) + alpha * runoff;
        dischargeData(time) = flow(end);
        flow = [0; slowStore - [slowStore(2:end); 0]];
        slowStore = slowStore - timeStep * diff(flow) + (1 - alpha) * runoff;
        dischargeData(time) = dischargeData(time) + flow(end);
        mainStoreData(time,:) = mainStore';
        fastStoreData(time,:) = fastStore';
        slowStoreData(time,:) = slowStore';
    end

    plot(1:(domainSize+1), [fastStoreData, zeros(nTimeSteps,1)]);
    plot(1:nTimeSteps, dischargeData);
end

