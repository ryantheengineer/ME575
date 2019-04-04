function [gradf] = gradient_f(x)
    gradf = [4*x(1)^3 - 4*x(1)*x(2) + 2*x(1) - 2;
         -2*x(1)^2 + 2*x(2)];
end