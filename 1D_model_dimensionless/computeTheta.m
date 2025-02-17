function [theta] = computeTheta(h, mod)
    theta = (h >= 0) .* mod.thetaS + (h < 0) .* ...
            (mod.thetaR + (mod.thetaS - mod.thetaR) .* ...
            (1 + (-mod.alpha .* h) .^ mod.n) .^ (-mod.m));
end
