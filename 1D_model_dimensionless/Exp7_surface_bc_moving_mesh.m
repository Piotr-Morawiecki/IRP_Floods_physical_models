mod.n = 1.282;
mod.m = 1 - 1/mod.n;
mod.alphaMvG = 3.367;
mod.ks = 5.4;
mod.z0 = 10;
mod.t0 = mod.z0 / mod.ks;
mod.alpha = mod.alphaMvG * mod.z0;
mod.thetaR = 0.115;
mod.thetaS = 0.388;

depth = 101;
dz = - repmat(1 / (depth - 1), depth - 1, 1);
z = cumsum([0; dz]);

totalTime = 0.0005;
nTimeSteps = 100;
dt = totalTime / nTimeSteps;
t = (dt * (0:(nTimeSteps-1)))';

maxIteration = 100;
maxError = 0.0001;

hData = zeros(depth, nTimeSteps);
zData = zeros(depth, nTimeSteps);
qData = zeros(depth + 1, nTimeSteps);
initialGroundwaterSurface = -1;
%initialGroundwaterSurface = 0;
hData(:, 1) = initialGroundwaterSurface - z;
zData(:, 1) = z;

flowIn = 1.5;
flowOut = 0;

mod.beta = 1;
%if flowIn > 0
%    mod.beta = flowIn / mod.ks;
%    mod.t0 = mod.z0 / (flowIn * mod.ks);
%else
%    mod.beta = 1;
%    mod.t0 = mod.z0 / mod.ks;
%end

b = 1;

surfaceReached = false;

for time = 2:nTimeSteps
    %disp(time);
    if ~surfaceReached && hData(1, time - 1) >= 0
        disp("Surface saturation reached after t = " + (time - 2) * dt);
        surfaceReached = true;
    end
    
    % Adapting mesh
    h0 = hData(:, time - 1);
    h = h0;
    gradH = [(h(1) - h(2)) / dz(1); ...
             (h(1:(end-2)) - h(3:end)) ./ (dz(1:(end-1)) + dz(2:end));...
             (h(end-1) - h(end)) / dz(end)];
    rho = sqrt(1 + b * gradH.^2);
    trapeziumArea = - (rho(1:(end - 1)) + rho(2:end)) / 2 .* dz;
    targetArea = sum(trapeziumArea) / (depth - 1);
    
    trapeziumId = 1;
    areaLeft = trapeziumArea(1);
    height = dz(1);
    base1 = rho(1);
    base2 = rho(2);
    for zId = 2:(depth - 1)
        currentArea = 0;
        while currentArea + areaLeft < targetArea
            currentArea = currentArea + areaLeft;
            trapeziumId = trapeziumId + 1;
            height = dz(trapeziumId);
            base1 = rho(trapeziumId);
            base2 = rho(trapeziumId + 1);
            areaLeft = trapeziumArea(trapeziumId);
        end
        area = targetArea - currentArea;
        gradient = (base2 - base1) / height;
        if gradient == 0
            deltaZ = area / base1;
        else
            deltaZ = ( - base1 + sqrt(base1 ^ 2 + 4 * area * gradient)) / ...
                     (2 * gradient);
        end
        base1 = base1 + gradient * deltaZ;
        height = height - deltaZ;
        z(zId) = z(trapeziumId + 1) + height / dz(trapeziumId) * ...
                 (z(trapeziumId) - z(trapeziumId + 1));
        h0(zId) = h(trapeziumId + 1) + height / dz(trapeziumId) * ...
                  (h(trapeziumId) - h(trapeziumId + 1));
        areaLeft = areaLeft - area;
    end
    dz = diff(z);
    middz = ( dz(1:(end-1)) + dz(2:end) ) / 2;
    
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
        upperRight = [diag(1 ./ dz); zeros(1, depth - 1)] - ...
                     [zeros(1, depth - 1); diag(1 ./ dz)];
        lowerLeft = [-diag(k ./ dz), zeros(depth - 1, 1)] + ...
                    [zeros(depth - 1, 1), diag(k ./ dz)];
        lowerRight = mod.beta * eye(depth - 1);
        matrixA = [upperLeft, upperRight ; lowerLeft, lowerRight];
        vectorB = [h0 .* sigma ./ dt; -k];
        vectorB(1) = vectorB(1) - flowIn / dz(1);
        vectorB(depth) = vectorB(depth) + flowOut / dz(end);

        solution = linsolve(matrixA, vectorB);
        error = norm((h - solution(1:depth)) / depth);
        iteration = iteration + 1;
        h = solution(1:depth);
    end
    if error > maxError
        warning("Warning: Max iterations exceeded with error " + error)
    end
    zData(:, time) = z;
    hData(:, time) = h;
    qData(:, time) = [-flowIn; solution((depth+1):end); -flowOut];
end

thetaData = computeTheta(hData, mod);
surfaceVolume = max(0, hData(1, :));
subsurfaceVolume = mean(thetaData, 1);
thetaData = effectiveSaturation(thetaData, mod);

figure
hold on
for time = 1:50
    plot3(zData(:, time), time * ones(depth, 1), hData(:, time));
end
hold off