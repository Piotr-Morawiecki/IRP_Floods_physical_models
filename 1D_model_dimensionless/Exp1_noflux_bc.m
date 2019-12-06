mod.n = 1.282;
mod.m = 1 - 1/mod.n;
mod.alphaMvG = 3.367;
mod.ks = 5.4;
mod.z0 = 10;
mod.t0 = mod.z0 / mod.ks;
mod.alpha = mod.alphaMvG * mod.z0;

depth = 20;
dz = - 1 / depth;
z = (dz * (1:depth))';

totalTime = 500;
nTimeSteps = 1000;
dt = totalTime / nTimeSteps;
t = (dt * (1:nTimeSteps))';

maxIteration = 10;
maxError = 0.000001;

hData = zeros(depth, nTimeSteps);
hData(:, 1) = -0.1;

for time = 2:nTimeSteps
    h0 = hData(:, time - 1);
    sigma = computeThetaDerivative(h0, mod);
    k = computeKr(h0, mod);
    k = (k(1:end-1) + k(2:end)) / 2;
    h = h0;
    error = 1;
    iteration = 1;
    while error > maxError && iteration < maxIteration
        upperLeft = diag(sigma ./ dt);
        upperRight = [eye(depth - 1); zeros(1, depth - 1)] ./ dz - ...
                     [zeros(1, depth - 1); eye(depth - 1)] ./ dz;
        lowerLeft = [-diag(k ./ dz), zeros(depth - 1, 1)] + ...
                    [zeros(depth - 1, 1), diag(k ./ dz)];
        lowerRight = eye(depth - 1);
        matrixA = [upperLeft, upperRight ; lowerLeft, lowerRight];
        vectorB = [h0 .* sigma ./ dt; -k];

        solution = linsolve(matrixA, vectorB);
        error = norm((h - solution(1:depth)) / depth);
        iteration = iteration + 1;
        h = solution(1:depth);
    end
    if error > maxError
        warning("Warning: Max iterations exceeded with error " + error)
    end
    hData(:, time) = h;
end

thetaData = computeTheta(hData, mod);
totalWaterVolume = mean(thetaData, 1);
figure;
subplot(3,1,1);
image(t, z, hData,'CDataMapping','scaled');
set(gca,'YDir','normal')
%image(hData,'CDataMapping','scaled');
title('phi');
colorbar;
subplot(3,1,2);
image(t, z, thetaData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('theta');
colorbar;
subplot(3,1,3);
title('total water volume');
plot(t, totalWaterVolume)