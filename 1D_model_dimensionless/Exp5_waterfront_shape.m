mod.n = 1.282;
mod.m = 1 - 1/mod.n;
mod.alphaMvG = 3.367;
mod.z0 = 10;
mod.alpha = mod.alphaMvG * mod.z0;
mod.thetaR = 0.115;
mod.thetaS = 0.388;

qG = 1.5;
c = qG / (mod.thetaS - mod.thetaR);

zspan = [0 5];
h0 = 0;
fun = @(z,h) - 1 / (computeThetaDerivative(h, mod) * c / computeKr(h, mod) + 1) ;
[z,h] = ode45(fun, zspan, h0);
figure;
subplot(2, 1, 1);
plot(z,h,'-')
ylabel("Hydraulic head");
xlabel("Distance before waterfront");
subplot(2, 1, 2);
plot(z,effectiveSaturation(computeTheta(h, mod), mod),'-')
ylabel("Effective saturation");
xlabel("Distance before waterfront");
