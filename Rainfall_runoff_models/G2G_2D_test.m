% size = dimensions of the gird
% gridSize = side length of each square [m]
% domain = shows if given square belongs to the domain
% river = shows if there is a river within given square

size = [10 10];

cat.dx = 1000;
cat.domain = true(size(1), size(2));
cat.river = false(size(1), size(2));
cat.river(7, 7:10) = true;
cat.gradient = zeros(size(1), size(2));
%cat.elevation = zeros(size(1), size(2));
%cat.conductivity = zeros(size(1), size(2));
%cat.specificYield = zeros(size(1), size(2));
%cat.vegetation = zeros(size(1), size(2));
cat.outlet = [7 10];

model.cl = 0.4;     % land surface speed [m/s]
model.cr = 0.5;     % river surface speed [m/s]
model.clb = 0.05;   % land sub-surface speed [m/s]
model.crb = 0.05;   % river sub-surface speed  [m/s]
model.rl = 0.005;   % land return flow factor
model.rr = 0.005;   % river return flow factor
model.cmax = 140;   % max store depth [mm]
model.kd = 0.00005; % drainage storage rate constant
model.beta = 3;

[surfaceAdjacency] = calculateRouting(cat, model, true);
[subsurfaceAdjacency] = calculateRouting(cat, model, false);

maxGradient = 10;
maxStorage = (1-cat.gradient./maxGradient)./model.cmax;

dt = 1;             % length of a time step [h]
nTimeSteps = 10;    % number of time steps during simulation

returnFlowFactor = (1 - cat.river) .* model.rl + cat.river .* model.rr;
thetaSurface = (1 - cat.river) .* model.cl .* dt ./ cat.dx + ...
               cat.river .* model.cr .* dt ./ cat.dx;
thetaSubSurface = (1 - cat.river) .* model.clb .* dt ./ cat.dx + ...
                  cat.river .* model.crb .* dt ./ cat.dx;

storage = zeros(size(1), size(2));
surfaceFlow = zeros(size(1), size(2));
subSurfaceFlow = zeros(size(1), size(2));
outletFlow = zeros(nTimeSteps, 2);

for time = 1:nTimeSteps
    evaporation = data.PE(time, :, :) .* ...
                  (1 - ((maxStorage - storage) ./ maxStorage) .^ 2);
    percipitation = data.percipitation(time, :, :);
    drainage = model.kd .* storage .^ beta;
    storage = storage + (percipitation - evaporation - drainage) .* dt;
    storage = max(storage, 0);
    
    directRunoff = max(0, store - maxStorage);
    storage = min(storage - maxStorage);
    
    surfaceRoute = (1 - thetaSurface) * surfaceRoute + ...
                   thetaSurface * (neighbours .* surfaceRoute + ...
                   directRunoff + returnFlowFactor .* subsurfaceRoute);
    subSurfaceRoute = (1 - thetaSubSurface) * surfaceRoute + ...
                    thetaSubSurface * (neighbours .* subSurfaceRoute + ...
                    directRunoff - returnFlowFactor .* subSurfaceRoute);
    outletFlow(time, :) = surfaceRoute(cat.outlet(1), cat.outlet(2)) + ...
                          thetaSubSurface(cat.outlet(1), cat.outlet(2));
end

outletTotalFlow = model.cr .* outletFlow(:, 1) + ...
                  model.crb .* outletFlow(:, 2); 