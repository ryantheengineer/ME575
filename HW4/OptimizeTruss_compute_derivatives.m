
    % ------------Starting point and bounds------------
    %design variables
    x0 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]; %starting point (all areas = 5 in^2)
    lb = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %lower bound
    ub = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20]; %upper bound
    global nfun;
    nfun = 0;
    
    % USE THIS FLAG TO CHOOSE THE DERIVATIVE TYPE:
    derivflag = 'forward';
%     derivflag = 'central';
    % derivflag = 'complex';

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
            [xopt, fopt, exitflag, output] = fmincon(@objfungrad_forward, x0, A, b, Aeq, beq, lb, ub, @confungrad_forward, options);
            xopt    %design variables at the minimum
            fopt    %objective function value at the minumum
            [f,~] = objfungrad_forward(xopt);
            [c,ceq,~,~] = confungrad_forward(xopt);
            c
            nfun
            
        case 'central'
            [xopt, fopt, exitflag, output] = fmincon(@objfungrad_central, x0, A, b, Aeq, beq, lb, ub, @confungrad_central, options);
            xopt    %design variables at the minimum
            fopt    %objective function value at the minumum
            [f,~] = objfungrad_central(xopt);
            [c,ceq,~,~] = confungrad_central(xopt);
            c
            nfun
            
        case 'complex'
            [xopt, fopt, exitflag, output] = fmincon(@objfungrad_complex, x0, A, b, Aeq, beq, lb, ub, @confungrad_complex, options);  
            
    end
    
   
    
    
    %% FORWARD DIFFERENCE FUNCTIONS %%
    % ------------Objective (forward difference)------------
    function [f,gradf] = objfungrad_forward(x)
        global nfun;
        
        %get data for truss from Data.m file
        Data;
        
        % insert areas (design variables) into correct matrix
        for i=1:nelem
            Elem(i,3) = x(i);
        end

        % call Truss to get weight and stresses
        [weight,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

        %objective function
        f = weight; %minimize weight
        
        % Gradient of the objective function
        if nargout > 1
            gradf = forward_partial(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
        end
        nfun = nfun + 1;
    end
        
    % ------------Non-linear Constraints (forward difference)------------
    function [c,ceq,DC,DCeq] = confungrad_forward(x)
        global nfun;
        Data;
        
        % insert areas (design variables) into correct matrix
        for i=1:nelem
            Elem(i,3) = x(i);
        end
    
        % call Truss to get stresses (constraints)
        [~,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

        c = zeros(10,1);
        for i = 1:10
            c(i) = stress(i);
        end

        % No nonlinear equality constraints
        ceq = [];

        % Gradient of the constraints
        if nargout > 2
            DC = forward_grad(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
            DCeq = [];
        end
        nfun = nfun + 1;
    end
    
    
    %% CENTRAL DIFFERENCE FUNCTIONS %%
    % ------------Objective (forward difference)------------
    function [f,gradf] = objfungrad_central(x)
        global nfun;
        
        %get data for truss from Data.m file
        Data;
        
        % insert areas (design variables) into correct matrix
        for i=1:nelem
            Elem(i,3) = x(i);
        end

        % call Truss to get weight and stresses
        [weight,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

        %objective function
        f = weight; %minimize weight
        
        % Gradient of the objective function
        if nargout > 1
            gradf = central_partial(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
        end
        nfun = nfun + 1;
    end
        
    % ------------Non-linear Constraints (forward difference)------------
    function [c,ceq,DC,DCeq] = confungrad_central(x)
        global nfun;
        Data;
        
        % insert areas (design variables) into correct matrix
        for i=1:nelem
            Elem(i,3) = x(i);
        end
    
        % call Truss to get stresses (constraints)
        [~,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

        c = zeros(10,1);
        for i = 1:10
            c(i) = stress(i);
        end

        % No nonlinear equality constraints
        ceq = [];

        % Gradient of the constraints
        if nargout > 2
            DC = central_grad(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
            DCeq = [];
        end
        nfun = nfun + 1;
    end