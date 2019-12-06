% simulation settings
domainSize = 10;
lengthStep = 1;
timeStep = 0.1;
nTimeSteps = 1000;

% callibrated parameters
alpha = 0.5;
fastFlowTime = 1;
slowFlowTime = 0.1;
maxStorage = repmat(1.1,domainSize,1);
rainfallRate = repmat(2,domainSize,1);

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