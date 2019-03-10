function [xopt, fopt, exitflag, output] = HW5externalopt()

    % ------------Starting point and bounds------------
    %var= Ps    Pf      delta  deltafactor    % Design variables
    x0 = [0.7,  0.001,  1,     0.9];      % Starting point
    ub = [0.99, 0.1,    2,     1.0];      % Upper bound
    lb = [0.5,  0.0001, 0.01,  0.1];      % Lower bound

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables
        Ps =  x(1);      % wire diameter
        Pf =  x(2);      % coil diameter
        delta =  x(3);      % number of coils
        deltafactor = x(4);      % unloaded height
        
        
        % Analysis functions
        pct_success = SimAnnealFunc(Ps,Pf,delta,deltafactor);
        
        
        % Objective function
        f = -pct_success;                           % Maximize the preload force
        
        % Inequality constraints
        c = zeros(1,1);
        c(1) = 0 - pct_success;
        
        % Equality constraints
        ceq = [];
        
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'iter-detailed');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    xopt
    fopt
    [~,c,~] = objcon(xopt);
    c
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end
end
