function [xopt, fopt, exitflag, output] = HW1()

    % ------------Starting point and bounds------------
    %var= d      D      n      hf           % Design variables
    x0 = [0.050, 0.500, 10.00, 1.500];      % Starting point
    ub = [0.250, 1.000, 100.0, 20.00];      % Upper bound
    lb = [0.001, 0.200, 1.000, 0.010];      % Lower bound

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables (THESE ARE CORRECT)
        d =  x(1);      % wire diameter
        D =  x(2);      % coil diameter
        n =  x(3);      % number of coils
        hf = x(4);      % unloaded height
        
        % Other analysis variables (THESE ARE CORRECT)
        h0 = 1.0;       % preload height
        delta0 = 0.4;   % deflection
        G = 12*10^6;    % psi
        Se = 45000;     % psi
        w = 0.18;
        Sf = 1.5;
        Q = 150000;     % psi
        
        
        % Analysis functions
        hdef = h0 - delta0;
        delta_max = (hf-h0) + delta0;
        delta_min = hf-h0;
        hs = n*d;   % solid height
        delta_hs = hf-hs;
        k = (G*d^4)/(8*(D^3)*n);    % spring stiffness
        F_preload = k*delta_min;    % preload force
        F_max = k*(hf-hdef);        % max working force
%         F_max = k*delta_hs;        % max working force
        F_s = k*delta_hs;          % force at the solid height
        K = (4*D-d)/(4*(D-d)) + 0.62*d/D;   % Wahl factor
        tau_max = (8*F_max*D*K)/(pi*d^3);   % stress due to F_max
        tau_min = (8*F_preload*D*K)/(pi*d^3);   % stress due to F_preload
%         tau_min = 0;
        tau_hs = (8*F_s*D*K)/(pi*d^3);
        tau_a = (tau_max - tau_min)/2;    % alternating stress
        tau_m = (tau_max + tau_min)/2;  % mean stress
        Sy = 0.44*Q/(d^w);    % yield strength in shear
        clash = hdef - hs;
        
        % Objective function
        f = -F_preload;
        
        % Inequality constraints
        c = zeros(9,1);
        c(1) = tau_hs - Sy;             % stress at hs <= Sy
        c(2) = D + d - 0.75;            % D + d <= 0.75
        c(3) = d - 0.2;                 % d <= 0.2
        c(4) = 0.01 - d;                % 0.01 <= d
        c(5) = 0.05 - clash;            % 0.05 <= clash
        c(6) = tau_a - (Se/Sf);         % tau_a <= Se/Sf
        c(7) = tau_a + tau_m -(Sy/Sf);  % tau_a + tau_m <= Sy/Sf
        c(8) = (D/d) - 16;
        c(9) = 4 - (D/d);
        
        Ddsum = D + d;   % Right now, the algorithm is taking this constraint right up to the boundary of 0.75
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