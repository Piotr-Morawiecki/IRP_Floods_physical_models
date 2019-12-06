function [t, z, hData, qData] = simulate1dFlow(beta, nSpatialSteps, ...
    totalTime, nTimeSteps, savingInterval)
load('DefaultModel.mat', 'mod');

depth = nSpatialSteps;
dz = - 1 / depth;
z = (dz * (1:depth))';

dt = totalTime / nTimeSteps;
t = savingInterval * dt * (0:(nTimeSteps / savingInterval - 1))';

maxIteration = 100;
maxError = 0.0001;

hData = zeros(depth, nTimeSteps/100);
qData = zeros(depth + 1, nTimeSteps/100);
initialGroundwaterSurface = -1;
hData(:, 1) = initialGroundwaterSurface - z;
h = hData(:, 1);

flowIn = beta;
flowOut = 0;

nextInfo = 100;
nextSave = max(2, savingInterval);

tic
for time = 2:nTimeSteps
    if time == nextInfo
        save('backup.mat', 'hData', 'qData');
        nextInfo = nextInfo + 100;
        timeMeasurement = toc;
        expectedTime = timeMeasurement * (nTimeSteps - time) / time;
        fprintf('Current time step: %d\n', time);
        fprintf('Progress made: %d%%\n', floor(100 * time / nTimeSteps));
        fprintf('Current time of computation: %f sec\n', timeMeasurement);
        fprintf('Estimated time of computation: ');
        if expectedTime > 3600
            fprintf('%f h\n', expectedTime / 3600);
        elseif expectedTime > 60
            fprintf('%f min\n', expectedTime / 60);
        else
            fprintf('%f s\n', expectedTime);
        end
        fprintf('----------\n');
    end
    h0 = h;
    error = 1;
    iteration = 1;
    while error > maxError && iteration < maxIteration
        sigma = computeThetaDerivative(h, mod);
        if h(1) >= 0
            sigma(1) = sigma(1) + 1;
        end
        k = computeKr(h, mod);
        k = (k(1:end-1) + k(2:end)) / 2;
        upperLeft = spdiags(sigma ./ dt, 0, depth, depth);
        upperRight = spdiags(repmat(1 / dz, depth, 1) * [1, -1], ...
                             [0 -1], depth, depth - 1);
        lowerLeft = spdiags(k / dz * [-1 1], [0 1], depth - 1, depth);
        lowerRight = spdiags(ones(depth - 1, 1), 0, ...
                             depth - 1, depth - 1);
        matrixA = [upperLeft, upperRight ; lowerLeft, lowerRight];
        vectorB = [h0 .* sigma ./ dt; -k];
        vectorB(1) = vectorB(1) - flowIn / dz;
        vectorB(depth) = vectorB(depth) + flowOut / dz;

        solution = mldivide(matrixA, vectorB);
        error = norm((h - solution(1:depth)) / depth);
        iteration = iteration + 1;
        h = solution(1:depth);
    end
    if error > maxError
        warning("Warning: Max iterations exceeded with error " + error);
    end
    if time == nextSave
        nextSave = nextSave + savingInterval;
        column = time / savingInterval;
        hData(:, column) = h;
        qData(:, column) = [-flowIn; solution((depth+1):end); -flowOut];
    end
end
end

