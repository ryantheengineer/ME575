%% Part a
% Check the given point for each inequality constraint
x = [-2.3723, -1.8364];

g2 = -x(1) - x(2)^2 + 1
g3 = -x(1) - x(2) + 1

%% Part c
clear;

fun = @prob3c;
x0 = [0,-0.75];
x = fsolve(fun,x0)