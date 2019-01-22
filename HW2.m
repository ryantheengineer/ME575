function [xopt, fopt, exitflag, output] = optimize_template()

    % ------------Starting point and bounds------------
    %var= V     D      d            % Design variables
    x0 = [];
    ub = [];
    lb = [];

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
        L = 15;             % length of pipeline, miles
        W = 12.67;          % flowrate of limestone, lbm/s
        a = 0.01;           % avg lump size of limestone before grinding, ft
        g = 32.17;          % acceleration due to gravity, ft/s^2
        rho_w = 62.4;       % density of water, lbm/ft^3
        gamma = 168.5;      % limestone density
        mu = 7.392*10^-4;   % viscosity of water lbm/(ft*s)
        
        
        % Analysis functions
        c = W/(Q_w + W);
        A = (pi/4)*D^2;
        rho = rho_w + c*(gamma-rho_w);
        mdot = rho*A*V;     % CHECK THIS TO MAKE SURE THE RIGHT RHO IS USED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Q = A*V;            % CHECK TO MAKE SURE THIS IS USED PROPERLY AS WELL!!!!!!!!!!!!!!!!!!!!!
        Pg = 218*W*((1/sqrt(d)) - (1/sqrt(a)));
        CdRpsq_calculated = 4*g*rho_w*(d^3)*((gamma-rho_w)/(3*mu^2));
        Cd = dragReynolds(CdRpsq_calculated);
        Rw = (rho_w*V*D)/mu;
        fw = fw_function(Rw);
        f = fw*((rho_w/rho) + 150*c*(rho_w/rho)*((g*D*(S-1))/((V^2)*sqrt(Cd)))^1.5);
        delta_p = (f*rho*L*V^2)/((D^2)*g_c);    % 
        
        
        % Objective function
        
        
        % Inequality constraints
        c = zeros(7,1);
        c(1) = ;            % 
        
        
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