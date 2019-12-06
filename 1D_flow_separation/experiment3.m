mod.thetaR = 0.115;
mod.thetaS = 0.388;
mod.n = 1.282;
mod.alphaMvG = 3.367;
mod.m = 1 - 1/mod.n;
mod.kSat = 5.4;
mod.grad = 0;

depth = 20;
dt = 1;
dz = -2;
maxIteration = 10;
maxError = 0.000001;

z = (dz * (1:depth))';
ks = computeKs(z, mod);
nTimeSteps = 1000;
hData = zeros(depth, nTimeSteps);
surfaceData = zeros(nTimeSteps, 1);
hData(:, 1) = -1;
rainfall = 0.002;
%surfaceData(1, 1) = 0;

for time = 2:nTimeSteps
    h0 = hData(:, time - 1);
    sigma = computeThetaDerivative(h0, mod);
    k = ks .* computeKr(h0, mod);
    
    drainage = -k(1) * (2 * (surfaceData(time - 1) - h0(1)) / dz + 1);
    surfaceData(time) = surfaceData(time - 1) + (rainfall - drainage) * dt;
    if surfaceData(time) < 0
        surfaceData(time) = 0;
        drainage = rainfall + surfaceData(time - 1) / dt;
    end
    
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
        vectorB(1) = vectorB(1) + drainage;

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
subsurfaceVolume = sum(thetaData, 1) * abs(dz);
totalWaterVolume = subsurfaceVolume + surfaceData';
figure;
subplot(3,1,1);
image(hData,'CDataMapping','scaled');
title('phi');
colorbar;
subplot(3,1,2);
image(thetaData,'CDataMapping','scaled');
title('theta');
colorbar;
subplot(3,1,3);
yyaxis left
plot(dt*(1:nTimeSteps), subsurfaceVolume)
yyaxis right
plot(dt*(1:nTimeSteps), surfaceData')
title('water volume per unit area');
legend('subsuface volume', 'surface volume');