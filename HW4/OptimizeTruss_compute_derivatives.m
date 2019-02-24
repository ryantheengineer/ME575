
    % ------------Starting point and bounds------------
    %design variables
    x0 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]; %starting point (all areas = 5 in^2)
    lb = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %lower bound
    ub = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20]; %upper bound
    global nfun;
    nfun = 0;

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % ------------Call fmincon------------
  
    options = optimoptions(@fmincon,'display','iter-detailed','Diagnostics','on','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
%     options = optimoptions(@fmincon,'display','iter-detailed','Diagnostics','on');
    [xopt, fopt, exitflag, output] = fmincon(@objfungrad, x0, A, b, Aeq, beq, lb, ub, @confungrad, options);  
%     [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);  
   
    xopt    %design variables at the minimum
    fopt    %objective function value at the minumum
%     [f, c, ceq] = objcon(xopt);
    [f,~] = objfungrad(xopt);
    [c,ceq,~,~] = confungrad(xopt);
    c
    nfun

    Data
    % ------------Objective------------
    function [f,gradf] = objfungrad(x)
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
    
        
        
    % ------------Non-linear Constraints------------
    function [c,ceq,DC,DCeq] = confungrad(x)
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
    
    
    
%     % ------------Objective and Non-linear Constraints------------
%     function [f, c, ceq] = objcon(x)
%         global nfun;
%         
%         %get data for truss from Data.m file
%         Data
%         
%         % insert areas (design variables) into correct matrix
%         for i=1:nelem
%             Elem(i,3) = x(i);
%         end
% 
%         % call Truss to get weight and stresses
%         [weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
% 
%         %objective function
%         f = weight; %minimize weight
%         
%         %inequality constraints (c<=0)
%         c = zeros(10,1);         % create column vector
%         for i=1:10
%             c(i) = abs(stress(i))-25000; % check stress both pos and neg         
%         end
%         
%         %equality constraints (ceq=0)
%         ceq = [];
%         nfun = nfun + 1;
% 
%     end
% 
%     % ------------Separate obj/con You may wish to change------------
%     function [f] = obj(x) 
%         [f, ~, ~] = objcon(x);
%     end
%     function [c, ceq] = con(x) 
%         [~, c, ceq] = objcon(x);
%     end