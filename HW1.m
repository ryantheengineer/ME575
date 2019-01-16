function [xopt, fopt, exitflag, output] = HW1()

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
        delta_max = (hf-h0) + delta0;
        delta_min = hf-h0;
        k = (G*d^4)/(8*D^3*n);  % spring stiffness
        F_preload = k*delta_min;   % preload force
        F_max = k*delta_max;       % max working force
        K = (4*D-d)/(4*(D-d)) + 0.62*d/D;   % Wahl factor
        tau_max = (8*F_max*D*K)/(pi*d^3);   % stress due to F_max
        tau_min = (8*F_preload*D*K)/(pi*d^3);   % stress due to F_preload
        tau_a = (tau_max - tau_min)/2;    % alternating stress
        tau_m = (tau_max + tau_min)/2;  % mean stress
        hs = n*d;   % solid height
        Sy = 0.44*Q/d^w;    % yield strength in shear
        
        
        % Objective function
        f = -F_preload;
        
        % Inequality constraints
        c = zeros(4,1); % HOW MANY INEQUALITY CONSTRAINTS?
        c(1) = stress_hs - Sy;  % stress at hs <= Sy
        c(2) = D + d - 0.75;
        c(3) = d - 0.2;
        c(4) = d + 0.01;    % not sure if this is the proper way to set this up
        
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