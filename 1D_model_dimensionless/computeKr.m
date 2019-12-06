function [kr] = computeKr(h, mod)
    theta = computeTheta(h, mod);
    Se = effectiveSaturation(theta, mod);
    kr = Se .^ 0.5 .* (1 - (1 - Se .^ (1 / mod.m)) .^ mod.m) .^ 2;
end