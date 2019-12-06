% simulation settings
domainSize = 10;
lengthStep = 1;
timeStep = 0.1;
nTimeSteps = 100;

% callibrated parameters
drainageTime = 0.1;
fastFlowTime = 1;
slowFlowTime = 0.01;
tresholdStorage = repmat(0.1, domainSize, 1);
maxStorage = repmat(1.1, domainSize, 1);
evaporationRate = repmat(2, domainSize, nTimeSteps);
% rainfallRate = repmat(2, domainSize, nTimeSteps);

avgRainfallRate = generateRainData(20,0.1,2,nTimeSteps);
rainfallRate = repmat(avgRainfallRate', domainSize, 1);

% types of stores
mainStore = zeros(domainSize,1);
fastStore = zeros(domainSize,1);
slowStore = zeros(domainSize,1);

% memory
mainStoreData = zeros(nTimeSteps, domainSize);
fastStoreData = zeros(nTimeSteps, domainSize);
slowStoreData = zeros(nTimeSteps, domainSize);
slowFlowData = zeros(nTimeSteps, domainSize+1);
fastFlowData = zeros(nTimeSteps, domainSize+1);
dischargeData = zeros(nTimeSteps, 1);

for time = 1:nTimeSteps
    drainage = (mainStore - tresholdStorage) * timeStep / drainageTime;
    drainage = timeStep * max(0, drainage);
    mainStore = mainStore - drainage + ...
        timeStep * (rainfallRate(:, time) - evaporationRate(:, time));
    directRunoff = max(0, mainStore - maxStorage);
    mainStore = max(0, min(maxStorage, mainStore));
    flow = [0; fastStore / fastFlowTime];
    fastFlowData(time,:) = flow;
    fastStore = fastStore - timeStep * diff(flow) + directRunoff;
    dischargeData(time) = flow(end);
    flow = [0; slowStore / slowFlowTime];
    slowFlowData(time,:) = flow;
    slowStore = slowStore - timeStep * diff(flow) + drainage;
    dischargeData(time) = dischargeData(time) + flow(end);
    mainStoreData(time,:) = mainStore';
    fastStoreData(time,:) = fastStore';
    slowStoreData(time,:) = slowStore';
end


hold on;
plot(1:nTimeSteps, -(avgRainfallRate*domainSize)/10, '-g');
plot(1:nTimeSteps, dischargeData, '-b');
legend('Rainfall rate','Discharge rate')