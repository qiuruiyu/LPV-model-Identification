function [y] = Add_noise(x, ratio)
    [n, ~] = size(x);
    e = randn(n, 1);
    v = filter([1], [1, -0.9], e);
    Px = var(x);
    Pv = var(v); 
    Pv_desired = ratio * Px;
    v = v * sqrt(Pv_desired/Pv);
    y = x + v; 
end

