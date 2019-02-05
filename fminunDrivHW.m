%---Example Driver program for fminun
    clear;

    global nobj ngrad
    nobj = 0; % counter for objective evaluations
    ngrad = 0.; % counter for gradient evaluations
    x0 = [1; 1; 1]; % starting point, set to be column vector
    algoflag = 1; % 1=steepest descent; 2=BFGS quasi-Newton
    stoptol = 1.e-3; % stopping tolerance, all gradient elements must be < stoptol  
    
    
    % ---------- call fminun----------------
    [xopt, fopt, exitflag] = fminun(@obj, @gradobj, x0, stoptol, algoflag);
   
    xopt
    fopt
    nobj
    ngrad
    
   
     % function to be minimized
     function [f] = obj(x)
        global nobj
        %test function 1
        f = 20 + 3*x(1) - 6*x(2) + 8*x(3) + 6*x(1)^2 - 2*x(1)*x(2) - ...
            x(1)*x(3) + x(2)^2 + 0.5*x(3)^2;
        nobj = nobj +1;
     end
     
    % get gradient as a column vector
     function [grad] = gradobj(x)
        global ngrad
        %gradient for test function 1
        grad(1,1) = 3 + 12*x(1) - 2*x(2) - x(3);
        grad(2,1) = -6 - 2*x(1) + 2*x(2);
        grad(3,1) = 8 - x(1) + x(3);
        ngrad = ngrad + 1;
     end
     
    