
    % ------------Starting point and bounds------------
    %design variables
    x0 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]; %starting point (all areas = 5 in^2)
    lb = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %lower bound
    ub = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20]; %upper bound
    global nfun;
    nfun = 0;

    % USE THIS FLAG TO CHOOSE THE DERIVATIVE TYPE:
%     derivflag = 'forward';
%     derivflag = 'central';
    derivflag = 'complex';
%     derivflag = 'default';

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

     % ------------Call fmincon------------
    options = optimoptions(@fmincon,'Algorithm','sqp','display','iter-detailed','Diagnostics','on');
    options = optimoptions(options,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

    switch derivflag
        case 'forward'
            [xopt, fopt, exitflag, output] = fmincon(@obj_forward, x0, A, b, Aeq, beq, lb, ub, @con_forward, options);
        case 'central'
            [xopt, fopt, exitflag, output] = fmincon(@obj_central, x0, A, b, Aeq, beq, lb, ub, @con_central, options);
        case 'complex'
            [xopt, fopt, exitflag, output] = fmincon(@obj_complex, x0, A, b, Aeq, beq, lb, ub, @con_complex, options);
        otherwise
            [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    end

    xopt    %design variables at the minimum
    fopt    %objective function value at the minumum
    [f, c, ceq] = objcon(xopt);
    c
    nfun

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        global nfun;

        %get data for truss from Data.m file
        Data;

        % insert areas (design variables) into correct matrix
        for i=1:nelem
            Elem(i,3) = x(i);
        end

        % call Truss to get weight and stresses
        [weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

        %objective function
        f = weight; %minimize weight

        %inequality constraints (c<=0)
        c = zeros(10,1);         % create column vector
        for i=1:10
            c(i) = sqrt((stress(i))^2)-25000; % check stress both pos and neg         
        end

        %equality constraints (ceq=0)
        ceq = [];
        nfun = nfun + 1;

    end

    % ------------Separate obj/con You may wish to change------------
    function [f] = obj(x) 
        [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x) 
        [~, c, ceq] = objcon(x);
    end

    
    %% Forward difference functions %%
    % ------------Forward objective------------
    function [f,gradf] = obj_forward(x)
        global nfun

        h = 0.001;
        n = numel(x);
        gradf = zeros(n,1);

        %objective function
        [f, ~, ~,] = objcon(x);

        % Gradient of the objective function
        if nargout > 1
            for i = 1:n
                xt = x;
                xt(i) = xt(i) + h;
                [fxt,~,~] = objcon(xt);
                gradf(i) = (fxt - f)/h;
            end
        end
        nfun = nfun + 1;
    end

    % ------------Forward constraints------------
    function [c,ceq,DC,DCeq] = con_forward(x)
        global nfun

        h = 0.001;
        n = numel(x);
        [~, c, ceq] = objcon(x);
        m = numel(c);

        DC = zeros(n,m);

        % Gradient of the constraints
        if nargout > 2
            % # of constraints
            for j = 1:m
                % # of values in each column
                for i = 1:n
                    xt = x;
                    xt(i) = xt(i) + h;
                    [~,cxt,~] = objcon(xt);
                    DC(i,j) = (cxt(j) - c(j))/h;
                end
                DCeq = [];
            end
        end
        nfun = nfun + 1;
    end

    
%% Central difference functions %%
    % ------------Central objective------------
    function [f,gradf] = obj_central(x)
        global nfun

        h = 0.001;
        n = numel(x);
        gradf = zeros(n,1);

        %objective function
        [f, ~, ~,] = objcon(x);

        % Gradient of the objective function
        if nargout > 1
            for i = 1:n
                xtplus = x;
                xtminus = x;
                xtplus(i) = xtplus(i) + h;
                xtminus(i) = xtminus(i) - h;
                [fxtplus,~,~] = objcon(xtplus);
                [fxtminus,~,~] = objcon(xtminus);
                gradf(i) = (fxtplus - fxtminus)/(2*h);
            end
        end
        nfun = nfun + 1;
    end

    % ------------Central constraints------------
    function [c,ceq,DC,DCeq] = con_central(x)
        global nfun

        h = 0.001;
        n = numel(x);
        [~, c, ceq] = objcon(x);
        m = numel(c);

        DC = zeros(n,m);

        % Gradient of the constraints
        if nargout > 2
            for j = 1:m
                for i = 1:n
                    xtplus = x;
                    xtminus = x;
                    xtplus(i) = xtplus(i) + h;
                    xtminus(i) = xtminus(i) - h;
                    [~,cxtplus,~] = objcon(xtplus);
                    [~,cxtminus,~] = objcon(xtminus);
                    DC(i,j) = (cxtplus(j) - cxtminus(j))/(2*h);
                end
                DCeq = [];
            end
        end
        nfun = nfun + 1;
    end


    %% Complex step functions %%
    % ------------Complex step objective------------
    function [f,grad_cs] = obj_complex(x)
        global nfun

        h = 0.0000001;
        n = numel(x);
        grad_cs = zeros(n,1);

        %objective function
        [f, ~, ~,] = objcon(x);

        % Gradient of the objective function
        if nargout > 1
            for i = 1:n
                z = complex(x,h);
                [fz,~,~] = objcon(z);
                grad_cs(i) = imag(fz)/h;
            end
        end
        nfun = nfun + 1;
    end

    % ------------Complex step constraints------------
    function [c,ceq,DC,DCeq] = con_complex(x)
        global nfun

        h = 0.0000001;
        n = numel(x);
        [~, c, ceq] = objcon(x);
        m = numel(c);

        DC = zeros(n,m);

        % Gradient of the constraints
        if nargout > 2
            for j = 1:m
                for i = 1:n
                    z = complex(x,h);
                    [~,cz,~] = objcon(z);
                    DC(i,j) = imag(cz(i))/h;
                end
                DCeq = [];
            end
        end
        nfun = nfun + 1;
    end