% simulation settings
domainSize = 10;
lengthStep = 1;
timeStep = 0.1;
nTimeSteps = 10000;

% callibrated parameters
drainageTime = 0.01;
fastFlowTime = 5;
slowFlowTime = 1000;
tresholdStorage = repmat(0.1, domainSize, 1);
maxStorage = repmat(1, domainSize, 1);
evaporationRate = repmat(2, domainSize, nTimeSteps);
% rainfallRate = repmat(2, domainSize, nTimeSteps);

avgRainfallRate = generateRainData(10,0.2,2,nTimeSteps);
rainfallRate = repmat(avgRainfallRate', domainSize, 1);


[fastFlowData, slowFlowData] = G2GLamb1999(domainSize, timeStep, ...
    nTimeSteps, drainageTime, fastFlowTime, slowFlowTime, ...
    tresholdStorage, maxStorage, rainfallRate, evaporationRate);

dischargeRate = fastFlowData + slowFlowData;

hold on;
plot(1:nTimeSteps, -(avgRainfallRate*domainSize)/10, '-k');
plot(1:nTimeSteps, slowFlowData, '-b', 'LineWidth', 2);
plot(1:nTimeSteps, fastFlowData + slowFlowData, '-r');
legend('Rainfall rate', 'Subsurface runoff', 'Surface runoff');