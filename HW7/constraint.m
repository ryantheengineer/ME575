function [g] = constraint(x)
    g = -(x(1) + 0.25)^2 + 0.75*x(2);
end