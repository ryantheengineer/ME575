% Homework 6, part 1 (KKT conditions)

% f = 4*x1 - 3*x2 + 2*x1^2 - 3*x1*x2 + 4*x2^2 (minimize)
% g1(x): 2*x1 - 1.5*x2 = 5

% KKT conditions
% (4 + 4*x1 - 3*x2) - lambda*2 = 0
% (-3 - 3*x1 + 8*x2) - lambda*(-1.5) = 0
% 2*x1 - 1.5*x2 - 5 = 0

% Matrix formulation of KKT conditions:
A = [4 -3   -2;
    -3  8   1.5;
     2 -1.5  0];

% x = [x1; x2; lambda];

b = [-4; 3; 5];

x = A\b
f = 4*x(1) - 3*x(2) + 2*x(1)^2 - 3*x(1)*x(2) + 4*x(2)^2

%% Part b
% If the constraint is changed...
b = [-4; 3; 5.1];

x = A\b
f = 4*x(1) - 3*x(2) + 2*x(1)^2 - 3*x(1)*x(2) + 4*x(2)^2
