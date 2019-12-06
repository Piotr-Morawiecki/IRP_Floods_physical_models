mod.n = 1.282;
mod.m = 1 - 1/mod.n;
mod.alphaMvG = 3.367;
mod.ks = 5.4;
mod.z0 = 10;
mod.t0 = mod.z0 / mod.ks;
mod.alpha = mod.alphaMvG * mod.z0;
mod.thetaR = 0.115;
mod.thetaS = 0.388;


dt = 0.0001;
timeMax = 0.2;
qSurf1 = runWavefront(1, dt, timeMax, mod);
qSurf2 = runWavefront(1.2, dt, timeMax, mod);
qSurf3 = runWavefront(1.5, dt, timeMax, mod);
qSurf4 = runWavefront(2, dt, timeMax, mod);
qSurf5 = runWavefront(3, dt, timeMax, mod);

figure;
   
subplot(1,2,1);
plot(qSurf1(:,1), qSurf1(:,3), ...
     qSurf2(:,1), qSurf2(:,3), ...
     qSurf3(:,1), qSurf3(:,3), ...
     qSurf4(:,1), qSurf4(:,3), ...
     qSurf5(:,1), qSurf5(:,3));
ylim([0,1.3]);
xlim([0,0.2]);
xlabel("Time")
ylabel("Groundwater flow")
legend("\beta = 1.01", "\beta = 1.2", "\beta = 1.5", "\beta = 2.0", ...
       "\beta = 3.0", 'Location', 'east');

subplot(1,2,2);
plot(qSurf1(:,1), qSurf1(:,2), ...
     qSurf2(:,1), qSurf2(:,2), ...
     qSurf3(:,1), qSurf3(:,2), ...
     qSurf4(:,1), qSurf4(:,2), ...
     qSurf5(:,1), qSurf5(:,2));
ylim([0,1.1]);
xlim([0,0.2]);
xlabel("Time")
ylabel("Ratio between surface and subsurface")
legend("\beta = 1.0", "\beta = 1.2", "\beta = 1.5", "\beta = 2.0", ...
       "\beta = 3.0", 'Location', 'east');

function [qSurf] = runWavefront(beta, dt, timeMax, mod)
zFront = 0;
time = 0;
qSurf = zeros(timeMax / dt, 3);
qSurf(:, 1) = linspace(0, timeMax, timeMax / dt);
while zFront > -1
    time = time + 1;
    dTheta = mod.thetaS - computeTheta(- 1 - zFront, mod);
    betaG = ((1 - dTheta) + sqrt((1 - dTheta)^2 + 4 * beta * dTheta)) / 2;
    qSurf(time, 3) = betaG;
    zFront = zFront - betaG / dTheta * dt;
end
qSurf(:, 2) = (beta - qSurf(:, 3)) ./ beta;
end