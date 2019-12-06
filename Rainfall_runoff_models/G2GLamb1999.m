function [fastFlowData, slowFlowData] = G2GLamb1999(domainSize, ...
    timeStep, nTimeSteps, drainageTime, fastFlowTime, slowFlowTime, ...
    tresholdStorage, maxStorage, rainfallRate, evaporationRate)

    % types of stores
    %mainStore = zeros(domainSize,1);
    mainStore = maxStorage;
    fastStore = zeros(domainSize,1);
    slowStore = zeros(domainSize,1);

    % memory
    %mainStoreData = zeros(nTimeSteps, domainSize);
    %fastStoreData = zeros(nTimeSteps, domainSize);
    %slowStoreData = zeros(nTimeSteps, domainSize);
    %slowFlowData = zeros(nTimeSteps, domainSize+1);
    %fastFlowData = zeros(nTimeSteps, domainSize+1);
    slowFlowData = zeros(nTimeSteps, 1);
    fastFlowData = zeros(nTimeSteps, 1);

    for time = 1:nTimeSteps
        drainage = (mainStore - tresholdStorage) * timeStep / drainageTime;
        drainage = timeStep * max(0, drainage);
        mainStore = mainStore - drainage + ...
            timeStep * (rainfallRate(:, time) - evaporationRate(:, time));
        directRunoff = max(0, mainStore - maxStorage);
        mainStore = max(0, min(maxStorage, mainStore));
        flow = [0; fastStore / fastFlowTime];
        %fastFlowData(time,:) = flow;
        fastStore = fastStore - timeStep * diff(flow) + directRunoff;
        fastFlowData(time) = flow(end);
        flow = [0; slowStore / slowFlowTime];
        %slowFlowData(time,:) = flow;
        slowStore = slowStore - timeStep * diff(flow) + drainage;
        slowFlowData(time) = flow(end);
        %mainStoreData(time,:) = mainStore';
        %fastStoreData(time,:) = fastStore';
        %slowStoreData(time,:) = slowStore';
    end
end