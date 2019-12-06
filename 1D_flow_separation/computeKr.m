function [kr] = computeKr(phi, mod)
    theta = computeTheta(phi, mod);
    Se = (theta - mod.thetaR) ./ (mod.thetaS - mod.thetaR);         % (8)
    kr = Se .^ 0.5 .* (1 - (1 - Se .^ (1 / mod.m)) .^ mod.m) .^ 2;  % (9)
end