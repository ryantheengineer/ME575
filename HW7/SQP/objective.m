function [f] = objective(x)
    f = x(1)^4 - 2*x(2)*x(1)^2 + x(2)^2 + x(1)^2 - 2*x(1) + 5;
end