clear;

model.cl = 0.4;     % land surface speed [m/s]
model.cr = 0.5;     % river surface speed [m/s]
model.clb = 0.05;   % land sub-surface speed [m/s]
model.crb = 0.05;   % river sub-surface speed  [m/s]
model.rl = 0.005;   % land return flow factor
model.rr = 0.005;   % river return flow factor
model.cmax = 140;   % max store depth [mm]
model.kd = 0.00005; % drainage storage rate constant
model.beta = 3;     % exponent in drainage equation 

% dim = dimensions of the gird
% dx = side length of each square [m]
% river = shows if there is a river within given square
% grad  = ...
% maxGrad = maximal gradient for ...
% outlet = ...
% surfaceAdjacency
% subSurfaceAdjacency

dim = [10 10];
cat.dx = 10;
cat.river = false(dim(1), dim(2));
cat.river(5, 1:5) = true;
cat.gradient = zeros(dim(1), dim(2));
cat.maxGrad = 10; 
cat.outlet = [5 1];
cat.surfaceAdjacency = calculateRouting(cat, model, true);
cat.subSurfaceAdjacency = calculateRouting(cat, model, false);
%cat.domain = true(dim(1), dim(2));
%cat.elevation = zeros(dim(1), dim(2));
%cat.conductivity = zeros(dim(1), dim(2));
%cat.specificYield = zeros(dim(1), dim(2));
%cat.vegetation = zeros(dim(1), dim(2));
saveas(image(cat.river, 'CDataMapping', 'scaled'), "catchment.png");

dt = 0.1;             % length of a time step [h]
nTimeSteps = 100000;  % number of time steps during simulation

data.potEvaporation = zeros(nTimeSteps, dim(1), dim(2));
%data.percipitation = 100*ones(nTimeSteps, dim(1), dim(2));
avgRainfallRate = generateRainData(100,0.1,4,nTimeSteps);
data.percipitation = repmat(avgRainfallRate, 1, dim(1), dim(2));

totalRainfall = avgRainfallRate * dim(1) * dim(2) * dt;
outletTotalFlow = runGrid2Grid(cat, model, data, nTimeSteps, dt);
              
plot(1:nTimeSteps, -totalRainfall/10, '-g');
hold on;
plot(1:nTimeSteps, outletTotalFlow, '-b');
legend('Rainfall rate','Discharge rate')
saveas(gcf,'outletFlow.png')