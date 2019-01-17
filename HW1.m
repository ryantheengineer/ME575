function [xopt, fopt, exitflag, output] = HW1()

    % ------------Starting point and bounds------------
    %var= d      D      n      hf           % Design variables
    x0 = [0.050, 0.500, 10.00, 1.500];      % Starting point
    ub = [0.250, 1.000, 100.0, 20.00];      % Upper bound
    lb = [0.010, 0.200, 1.000, 0.010];      % Lower bound
    ub = [0.200, 1.000, 100.0, 20.00];      % Upper bound
    lb = [0.01, 0.200, 1.000, 0.010];      % Lower bound

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables (THESE ARE OK)
        d =  x(1);      % wire diameter
        D =  x(2);      % coil diameter
        n =  x(3);      % number of coils
        hf = x(4);      % unloaded height
        
        % Other analysis variables (THESE ARE OK)
        h0 = 1.0;               % preload height, in
        delta0 = 0.4;           % deflection, in
        hdef = h0 - delta0;     % deflection height, in
        G = 12*10^6;            % psi
        Se = 45000;             % psi
        w = 0.18;
        Sf = 1.5;
        Q = 150000;             % psi
        delta_f = 0;            % in
        
        % Analysis functions
        k = (G*d^4)/(8*(D^3)*n);            % spring stiffness
        delta_p = hf-h0;                    % deflection at preload
        delta_def = delta0 + (hf-h0);       % greatest working deflection
        hs = n*d;                           % solid height
        delta_s = hf - hs;                  % deflection at solid height
        F_f = k*delta_f;                    % full height force (zero)
        F_p = k*delta_p;                    % preload force
        F_def = k*delta_def;                % force at full deflection
        F_s = k*delta_s;                    % force at solid height
        K = ((4*D)-d)/(4*(D-d))+0.62*(d/D); % Wahl factor
        tau_f = ((8*F_f*D)/(pi*d^3))*K;     % stress at full height (zero)
        tau_p = ((8*F_p*D)/(pi*d^3))*K;     % stress at preload height
        tau_def = ((8*F_def*D)/(pi*d^3))*K; % stress at full deflection
        tau_s = ((8*F_s*D)/(pi*d^3))*K;     % stress at solid height
        tau_max = tau_def;                  % max stress (assumed at max deflection)
        tau_min = tau_p;                    % min stress (assumed at preload deflection)
        tau_m = (tau_max + tau_min)/2;      % mean stress
        tau_a = (tau_max - tau_min)/2;      % alternating stress
        Sy = 0.44*(Q/(d^w));                % yield stress
        dratio = D/d;                       % ratio of diameters
        dsum = D + d;                       % sum of diameters
        clash = hdef - hs;                  % clash allowance
        Sefratio = Se/Sf;                   % endurance limit to safety factor
        Syfratio = Sy/Sf;                   % yield strength to safety factor
        tau_amsum = tau_a + tau_m;          % sum of alternating and mean stresses
        
        % Objective function
        f = -F_p;                           % Maximize the preload force
        
        % Inequality constraints
        c = zeros(7,1);
        c(1) = tau_a - Sefratio;            % tau_a <= Se/Sf
        c(2) = tau_amsum - Syfratio;        % tau_a + tau_m <= Sy/Sf
        c(3) = dratio - 16;                 % D/d <= 16
        c(4) = 4 - dratio;                  % 4 <= D/d
        c(5) = dsum - 0.75;                 % D + d <= 0.75
        c(6) = 0.05 - clash;                % 0.05 <= clash allowance
        c(7) = tau_s - Sy;                  % tau_s <= Sy
        
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