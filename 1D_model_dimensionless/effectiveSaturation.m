function [se] = effectiveSaturation(theta, mod)
    se = (theta - mod.thetaR) ./ (mod.thetaS - mod.thetaR);
end

