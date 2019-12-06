function [dTheta] = computeThetaDerivative(phi, mod)
    dTheta = (phi<0) .* mod.alphaMvG .* (mod.thetaS - mod.thetaR) .* ...
             mod.n .* (-mod.alphaMvG .* phi) .^ (mod.n - 1) .* mod.m .* ...
             ((-mod.alphaMvG .* phi) .^ mod.n + 1) .^ (-mod.m - 1);
end