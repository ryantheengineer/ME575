% HW7_SQP.m: Perform the next two iterations of the SQP example problem in
% the class notes, picking up where the example leaves off.

% Min f(x) = x1^4 - 2*x2*x1^2 + x2^2 + x1^2 - 2*x1 + 5
% s.t. g(x) = -(x1 + 0.25)^2 + 0.75*x2 >= 0
% Starting from the point [-1,4]

% (For the homework problem, this actually means starting from the point:
% transpose(x) = [0.2533,-0.1583]


% Gradient of the objective
x = [-1; 4];
f = objective(x)
gradf = gradient_f(x)
g = constraint(x)
gradg = gradient_g(x)
HessL = eye(2)




%% Functions
function [f] = objective(x)
    f = x(1)^4 - 2*x(2)*x(1)^2 + x(2)^2 + x(1)^2 - 2*x(1) + 5;
end

function [gradf] = gradient_f(x)
    gradf = [4*x(1)^3 - 4*x(1)*x(2) + 2*x(1) - 2;
         -2*x(1)^2 + 2*x(2)];
end

function [g] = constraint(x)
    g = -(x(1) + 0.25)^2 + 0.75*x(2);
end

function [gradg] = gradient_g(x)
    gradg = [-2*x(1) - 0.5;
             0.75];
end

% function [fa] = 