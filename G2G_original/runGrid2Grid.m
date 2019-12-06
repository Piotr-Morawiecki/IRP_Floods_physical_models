function [outletTotalFlow] = runGrid2Grid(cat, model, data, nTimeSteps, dt)
    dim = size(cat.river);
    area = dim(1) * dim(2);

    maxStorage = (1-cat.gradient./cat.maxGrad).*model.cmax;

    returnFlowFactor = (1 - cat.river) .* model.rl + cat.river .* model.rr;
    thetaSurface = (1 - cat.river) .* model.cl .* dt ./ cat.dx + ...
                   cat.river .* model.cr .* dt ./ cat.dx;
    thetaSubSurface = (1 - cat.river) .* model.clb .* dt ./ cat.dx + ...
                      cat.river .* model.crb .* dt ./ cat.dx;

    storage = zeros(dim(1), dim(2));
    surfaceFlow = zeros(dim(1), dim(2));
    subSurfaceFlow = zeros(dim(1), dim(2));
    outFlow = zeros(nTimeSteps, 2);
    storageMemory = zeros(nTimeSteps, dim(1), dim(2));
    %summary = zeros(nTimeSteps, 3);

    for time = 1:nTimeSteps
        evaporation = squeeze(data.potEvaporation(time, :, :)) .* ...
                      (1 - ((maxStorage - storage) ./ maxStorage) .^ 2);
        percipitation = squeeze(data.percipitation(time, :, :));
        drainage = model.kd .* storage .^ model.beta;
        storageMemory(time, :, :) = storage;
        storage = storage + (percipitation - evaporation - drainage) .* dt;
        storage = max(storage, 0);

        directRunoff = max(0, storage - maxStorage);
        storage = storage - directRunoff;

        inflow = reshape(cat.surfaceAdjacency * ...
                         reshape(surfaceFlow, area, 1), dim(1), dim(2));
        surfaceFlow = (1 - thetaSurface) .* surfaceFlow + ...
                      thetaSurface .* inflow + directRunoff + ...
                      returnFlowFactor .* subSurfaceFlow * dt;
        inflow = reshape(cat.subSurfaceAdjacency * ...
                         reshape(subSurfaceFlow, area, 1), dim(1), dim(2));
        subSurfaceFlow = (1 - thetaSubSurface) .* subSurfaceFlow + ...
                         thetaSubSurface .* inflow + drainage * dt - ...
                         returnFlowFactor .* subSurfaceFlow * dt;
        outFlow(time, :) = [surfaceFlow(cat.outlet(1), cat.outlet(2)) , ...
                           subSurfaceFlow(cat.outlet(1), cat.outlet(2))];
        %summary(time, 1) = sum(storage, 'all');
        %summary(time, 2) = sum(surfaceFlow, 'all');
        %summary(time, 3) = sum(subSurfaceFlow, 'all');
    end

    outletTotalFlow = thetaSurface(cat.outlet(1), cat.outlet(2)) .* ...
                      outFlow(:, 1) + thetaSubSurface(cat.outlet(1), ...
                      cat.outlet(2)) .* outFlow(:, 2);
end

