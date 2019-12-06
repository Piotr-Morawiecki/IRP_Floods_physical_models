function [phi] = theta2phi(theta, mod)
    if theta == mod.thetaS
        phi = 0;
    else
        phi = -1 / mod.alphaMvG * (((theta - mod.thetaR) / ...
            (mod.thetaS - mod.thetaR)) ^ (-1 / mod.m) - 1) ^ (1 / mod.n);
    end
end

