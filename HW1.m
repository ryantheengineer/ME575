function [xopt, fopt, exitflag, output] = optimize_template()

    % ------------Starting point and bounds------------
    %var= d     D    n    hf    %design variables
    x0 = [0.05, 0.5, 10,  1.5];
    ub = [0.2,  0.74,100, 20];
    lb = [0.01, 0.55,1,   0.01];

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables
        d = x(1);   % wire diameter
        D = x(2);   % coil diameter
        n = x(3);   % number of coils
        hf = x(4);  % 
        
        % Other analysis variables
        h0 = 1.0;       % preload height
        delta0 = 0.4;   % deflection
        G = 12*10^6;    % psi
        Se = 45000;     % psi
        w = 0.18;
        Sf = 1.5;
        Q - 150000;     % psi
        
        % Analysis functions
        
        % Objective function
        
        % Inequality constraints
        c = zeros(3,1);
        c(1) = stress_hs - Sy;  % stress at hs <= Sy
        c(2) = 
        
        % Equality constraints
        ceq = [];
        
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'iter-detailed');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end
end