
beta = 0.1;
nSpatialSteps = 1000;
totalTime = 5.04;
nTimeSteps = 10000;
savingInterval = 10;

%{
beta = 1.2;
nSpatialSteps = 30*1000;
totalTime = 0.45;
nTimeSteps = 30*10000;
savingInterval = 100;
%}

tic
[t, z, hData, qData] = simulate1dFlow(beta, nSpatialSteps, totalTime, ...
                                      nTimeSteps, savingInterval);
toc

save('result4.mat')

%hData = hData(:, 1:7976);
%t = t(1:7976);
plotSimulationResults(t, z, hData, qData);
plotRegionsForHighRainfall(z, t, hData, 10e-5);


%region1H = hData(1, 9000);
%epsilon = 10^(-3);
%plotRegionBoundaries(z, t, hData, region1H, epsilon)
