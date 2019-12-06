% https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2015WR018508
% Parameter values were taken from Table 1.
mod.thetaR = 0.115;
mod.thetaS = 0.388;
mod.n = 1.282;
mod.alphaMvG = 3.367;
mod.m = 1 - 1/mod.n;
mod.kSat = 5.4;
mod.grad = 0;

%h = -0.001*2.^(0:20);
%theta = arrayfun(@(x) computeTheta(x, mod), h);
%kr = arrayfun(@(x) computeKr(x, mod), theta);
%semilogx(h, kr, h, theta);

maxDepth = 2;
nTimeSteps = 2;
phi = zeros(nTimeSteps, maxDepth);
phi(1, 1:maxDepth) = -1;
dt = 0.1;
dz = 1;
z = -dz * (1:maxDepth);
ks = arrayfun(@(x) computeKs(x, mod), z);

for time = 2:nTimeSteps
    phiGrad = (phi(time-1, 1:(maxDepth-1)) - phi(time-1, 2:maxDepth)) / dz;
    %thetaDeriv = arrayfun(@(x) computeThetaDerivative(x, mod), ...
    %                      phi(time - 1, :));
    theta = arrayfun(@(x) computeTheta(x, mod), phi(time - 1, :));
    thetaDeriv = diff(theta) ./ dz;
    k = ks .* arrayfun(@(x) computeKr(computeTheta(x, mod), mod), ...
                       phi(time - 1, :));
    k = (k(1:(maxDepth-1)) + k(2:maxDepth)) / 2;
    v = [0, -k .* (diff(phi(time - 1, :)) ./ dz + 1), 0];
    phi(time, :) = phi(time - 1, :) + diff(v) .* dt ./ dz ./ thetaDeriv;
    %phi(time, 2:maxDepth) = phi(time - 1, 2:maxDepth) + ...
    %                        ( thetaDerivative .* k .* (phiGrad + 1)) .* dt;
end

thetaArray = arrayfun(@(x) computeTheta(x, mod), phi');
totalWaterVolume = mean(thetaArray, 1);
figure;
subplot(3,1,1);
image(phi','CDataMapping','scaled')
title('phi')
subplot(3,1,2);
image(thetaArray,'CDataMapping','scaled')
title('theta')
subplot(3,1,3);
plot(totalWaterVolume)
