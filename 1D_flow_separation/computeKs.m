function [ks] = computeKs(z, mod)
    ks = mod.kSat + mod.grad .* z;
end