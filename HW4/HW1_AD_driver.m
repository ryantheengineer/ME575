% This routine calls the analysis program HW1_AD.m, which automatically
% computes derivatives

% Set initial values of variables
x1 = 10;       % tau_a
x2 = 100;       % tau_amsum
x3 = 2;         % dratio
x4 = 1;         % dsum
x5 = 0.03;      % clash
x6 = 100;       % tau_s

% Call the analysis routine
[Functions, Jacobian] = HW1_AD(x1,x2,x3,x4,x5,x6);

% Print values
Functions
Jacobian
