function [kr] = computeKr2(h, mod)
    ah = - mod.alpha .* h;
    kr = (h < 0) .* (1 - ah .^ (mod.n - 1) .* (1 + ah .^ mod.n) .^ ...
         (-mod.m)) .^ 2 ./ (1 + ah .^ mod.n) .^ (mod.m / 2) + (h >= 0);
end