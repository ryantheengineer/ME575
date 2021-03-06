function [xopt, fopt, exitflag, output] = HW2()

    % ------------Starting point and bounds------------
    %var= V     D      d            % Design variables
    x0 = [14.96,0.4,   0.005];
    ub = [30,   0.5,   0.007];
    lb = [0.01, 0.1,   0.0005];

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        % Design variables
        V =  x(1);      % average flow velocity, ft/sec
        D =  x(2);      % internal pipe diameter, ft
        d =  x(3);      % avg limestone particle size after grinding, ft

        % Other analysis variables
        L = 15*5280;        % length of pipeline, feet
        W = 12.67;          % flowrate of limestone, lbm/s
        a = 0.01;           % avg lump size of limestone before grinding, ft
        g = 32.17;          % acceleration due to gravity, ft/s^2
        rho_w = 62.4;       % density of water, lbm/ft^3
        gamma = 168.5;      % limestone density
        S = gamma/rho_w;    % limestone specific gravity
        mu = 7.392*10^-4;   % viscosity of water lbm/(ft*s)
        gc = 32.17;         % conversion factor between lbf and lbm
        
        % Analysis functions
        C = 4*W/(pi*gamma*V*(D^2))
        Area = (pi/4)*D^2
        rho = rho_w + C*(gamma-rho_w)
        Pg = (218*W*((1/sqrt(d)) - (1/sqrt(a))))/550
        CdRpsq_calculated = 4*g*rho_w*(d^3)*((gamma-rho_w)/(3*mu^2))
        Cd = dragReynolds(CdRpsq_calculated)
        Rw = (rho_w*V*D)/mu
        fw = fw_function(Rw)
        F = fw*((rho_w/rho) + 150*C*(rho_w/rho)*((g*D*(S-1))/((V^2)*sqrt(Cd)))^1.5)
        delta_p = (F*rho*L*V^2)/(2*D*gc) 
        Qslurry = Area*V
        Pf = (delta_p*Qslurry)/550     % pump power, hp
        Vc = ((40*g*C*(S-1)*D)/sqrt(Cd))^0.5
        horsepower = Pf + Pg
        cost = NPV(horsepower,0.07,7) % develop cost function here
        
        % Objective function
        f = cost;               % Minimize the total cost
        
        % Inequality constraints
        c = zeros(2,1);
        c(1) = (1.1*Vc) - V;        % 1.1*Vc <= V
        c(2) = C - 0.4;             % C <= 0.4        
        
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