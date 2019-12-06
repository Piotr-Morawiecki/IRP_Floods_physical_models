mod.n = 1.282;
mod.m = 1 - 1/mod.n;
mod.alphaMvG = 3.367;
mod.ks = 5.4;
mod.z0 = 10;
mod.t0 = mod.z0 / mod.ks;
mod.alpha = mod.alphaMvG * mod.z0;
mod.thetaR = 0.115;
mod.thetaS = 0.388;

refinedZone = 0.02;
refinedZoneDepth = 400;
sparseZoneDepth = 20;
depth = refinedZoneDepth + sparseZoneDepth;
dz = - [repmat(refinedZone / refinedZoneDepth, refinedZoneDepth, 1); ...
        repmat((1 - refinedZone) / sparseZoneDepth, sparseZoneDepth, 1)];
middz = ( dz(1:(end-1)) + dz(2:end) ) / 2;
z = cumsum(dz);

totalTime = 0.005;
nTimeSteps = 2000;
dt = totalTime / nTimeSteps;
t = (dt * (0:(nTimeSteps-1)))';

maxIteration = 100;
maxError = 0.0001;

hData = zeros(depth, nTimeSteps);
qData = zeros(depth + 1, nTimeSteps);
initialGroundwaterSurface = -1;
%initialGroundwaterSurface = 0;
hData(:, 1) = initialGroundwaterSurface - z;

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

surfaceReached = false;

tic
for time = 2:nTimeSteps
    %disp(time);
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
        %tic
        upperLeft = spdiags(sigma ./ dt, 0, depth, depth);
        upperRight = spdiags([1 ./ dz(1:(end - 1)), -1 ./ dz(2:end)], ...
                             [0 -1], depth, depth - 1);
        lowerLeft = spdiags(k ./ middz * [-1 1], [0 1], depth - 1, depth);
        lowerRight = spdiags(repmat(mod.beta, depth - 1, 1), 0, ...
                             depth - 1, depth - 1);
        %matrixA = sparse([upperLeft, upperRight ; lowerLeft, lowerRight]);
        %toc
        matrixA = [upperLeft, upperRight ; lowerLeft, lowerRight];
        vectorB = [h0 .* sigma ./ dt; -k];
        vectorB(1) = vectorB(1) - flowIn / dz(1);
        vectorB(depth) = vectorB(depth) + flowOut / dz(end);
        
        %tic
        solution = mldivide(matrixA, vectorB);
        %toc
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
toc

thetaData = computeTheta(hData, mod);
surfaceVolume = max(0, hData(1, :));
subsurfaceVolume = mean(thetaData, 1);
thetaData = effectiveSaturation(thetaData, mod);

%figure
%hold on
%for time = (1:100)*100
%    plot3(z, time * ones(size(z)), hData(:, time))
%end
%hold off

figure;
hold on
region2_top = z(find(hData(1:(end - 1), end) > 0, 1, 'last'));
region2_bottom = z(hData(:, end) == min(hData(:, end)));
plot([-refinedZone 0], [refinedZone-1 -1], 'b');
plot([-refinedZone 0], [- h(1,end) / region2_top * (- refinedZone) + ...
     h(1, end), hData(1, end)], 'c');
plot(z(1:refinedZoneDepth), hData(1:refinedZoneDepth, end), 'r');
plot([region2_top region2_top], [-1 0.1], 'g');
plot([region2_bottom region2_bottom], [-1 0.1], 'g');
axis([-refinedZone 0 -1 0.1])
xlabel('Depth (z)')
ylabel('Pressure head (h)')
legend('Region I solution', 'Region III solution', 'Full solution', ...
       'Region boundaries');
hold off

%plot(h(1:refinedZoneDepth), hData(1:refinedZoneDepth, end) - (-h(1:refinedZoneDepth)-1))