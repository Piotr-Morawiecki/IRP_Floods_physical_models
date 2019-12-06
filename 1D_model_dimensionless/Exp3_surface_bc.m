mod.n = 1.282;
mod.m = 1 - 1/mod.n;
mod.alphaMvG = 3.367;
mod.ks = 5.4;
mod.z0 = 10;
mod.t0 = mod.z0 / mod.ks;
mod.alpha = mod.alphaMvG * mod.z0;
mod.thetaR = 0.115;
mod.thetaS = 0.388;

depth = 100;
dz = - 1 / depth;
z = (dz * (1:depth))';

totalTime = 0.1;
nTimeSteps = 10000;
dt = totalTime / nTimeSteps;
t = (dt * (0:(nTimeSteps-1)))';

maxIteration = 100;
maxError = 0.0001;

hData = zeros(depth, nTimeSteps);
qData = zeros(depth + 1, nTimeSteps);
initialGroundwaterSurface = -1;
%initialGroundwaterSurface = 0;
hData(:, 1) = initialGroundwaterSurface - z;

flowIn = 1.1;
flowOut = 0;

mod.beta = 1;
%if flowIn > 0
%    mod.beta = flowIn / mod.ks;
%    mod.t0 = mod.z0 / (flowIn * mod.ks);
%else
%    mod.beta = 1;
%    mod.t0 = mod.z0 / mod.ks;
%end

surfaceReached = false;

for time = 2:nTimeSteps
    if ~surfaceReached && hData(1, time - 1) >= 0
        disp("Surface saturation reached after t = " + (time - 2) * dt);
        surfaceReached = true;
    end
    h0 = hData(:, time - 1);
    h = h0;
    error = 1;
    iteration = 1;
    while error > maxError && iteration < maxIteration
        sigma = computeThetaDerivative(h, mod);
        if h(1) >= 0
            sigma(1) = sigma(1) + 1;
        end
        k = computeKr(h, mod);
        k = (k(1:end-1) + k(2:end)) / 2;
        upperLeft = diag(sigma ./ dt);
        upperRight = [eye(depth - 1); zeros(1, depth - 1)] ./ dz - ...
                     [zeros(1, depth - 1); eye(depth - 1)] ./ dz;
        lowerLeft = [-diag(k ./ dz), zeros(depth - 1, 1)] + ...
                    [zeros(depth - 1, 1), diag(k ./ dz)];
        lowerRight = mod.beta * eye(depth - 1);
        matrixA = [upperLeft, upperRight ; lowerLeft, lowerRight];
        vectorB = [h0 .* sigma ./ dt; -k];
        vectorB(1) = vectorB(1) - flowIn / dz;
        vectorB(depth) = vectorB(depth) + flowOut / dz;

        solution = linsolve(matrixA, vectorB);
        error = norm((h - solution(1:depth)) / depth);
        iteration = iteration + 1;
        h = solution(1:depth);
    end
    if error > maxError
        warning("Warning: Max iterations exceeded with error " + error)
    end
    hData(:, time) = h;
    qData(:, time) = [-flowIn; solution((depth+1):end); -flowOut];
end

thetaData = computeTheta(hData, mod);
surfaceVolume = max(0, hData(1, :));
subsurfaceVolume = mean(thetaData, 1);
thetaData = effectiveSaturation(thetaData, mod);
figure;

subplot(4,1,1);
image(t, z, hData,'CDataMapping','scaled');
set(gca,'YDir','normal')
%image(hData,'CDataMapping','scaled');
title('Hydraulic head');
colorbar;

subplot(4,1,2);
image(t, z, thetaData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('Effective saturation');
colorbar;

subplot(4,1,3);
image(t, [0;z], qData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('Flow');
colorbar;

subplot(4,1,4);

yyaxis left
plot(t, subsurfaceVolume)
yyaxis right
plot(t, surfaceVolume)
title('water volume per unit area');
legend('subsuface volume', 'surface volume');
title('Total water volume');