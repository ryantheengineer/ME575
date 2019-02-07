
     function [xopt, fopt, exitflag] = fminun(obj, gradobj, x0, stoptol, algoflag)
        global nobj
        % get function and gradient at starting point
        [n,~] = size(x0); % get number of variables
        f = obj(x0)
        grad = gradobj(x0)
        x = x0;       
     
        if (algoflag == 1)      % steepest descent
            while nobj < 20000
                % set an approach tolerance for Newton step
                apprchtol = 1;

                %set starting step length
                alpha = 0.12;
                s = srchsd(grad);
                alpha_star = minimizing_step(obj,s,x,f,alpha);
                [xnew,fnew] = take_step(obj,x,alpha_star,s);
                % Check to see if the gradient is within the desired tolerance
                grad = gradobj(xnew)
                if abs(grad(1))>apprchtol || abs(grad(2))>apprchtol || abs(grad(3))>apprchtol % if the gradient is not smaller than stoptol
                    x = xnew;
                    f = fnew;
                    continue
                else
                    x = xnew;
                    f = fnew;
                    [xnew,fnew] = Newton_quad(obj,gradobj,x);
                    exitflag = 1;
                    break
                end
%                 if abs(grad(1))>apprchtol || abs(grad(2))>apprchtol || abs(grad(3))>apprchtol % if the gradient is not smaller than stoptol
%                     x = xnew;
%                     f = fnew;
%                     continue
%                 else
%                     x = xnew;
%                     f = fnew;
%                     [xnew,fnew] = Newton_quad(obj,gradobj,x);
%                     exitflag = 1;
%                     break
%                 end
            end
            
        end
        
%         if (algoflag == 2)      % BFGS quasi-Newton
%             s = srchbfgs(grad)
%         end
        
%             if all(abs(grad)>stoptol)   % if the gradient is not smaller than the

        % while the number of function calls is less than our set maximum
%         while nobj < 1000   
            
            
            
                
        grad = gradobj(xnew)
        if grad > stoptol
            error('Not under stoptol. Adjust apprchtol and run again.');
        end
        xopt = xnew;
        fopt = fnew;
        
        exitflag = 0;   % This is automatically set to stop the solver from iterating for a long time, and it's why xopt doesn't result in a minimum right now
     end
     
     % get steepest descent search direction as a column vector
     function [s] = srchsd(grad)
        mag = sqrt(grad'*grad);
        s = -grad/mag;
     end
     
     function [alpha_star] = minimizing_step(obj,s,x,f,alpha)
        % Create holding vectors for the f and alpha values
        f_alpha = [f,alpha];
        count = 0;
        
        while count < 1000   % this might not be the right limit
            % take a step
            xnew = x + alpha*s;
            fnew = obj(xnew);
            
            % check if the value of fnew is less than f (is the function
            % still decreasing?):
            if fnew < f
                x = xnew;   % reset x for the next step
                alpha = 2*alpha;    % double alpha for the next step
                f_alpha_temp = [fnew,alpha];
                f_alpha = [f_alpha; f_alpha_temp];
                f = fnew;   % set the new f for a point of comparison on the next time through the loop
                count = count + 1;
                continue
            else
                % if the function begins to increase, then curve fit the
                % data with a parabola, and step to the minimum of the
                % parabola.
                alpha = alpha/2;
                f_alpha_temp = [fnew,alpha];
                f_alpha = [f_alpha; f_alpha_temp];
                  
                break
            end
        end
        
        % Reorder f_alpha so that the "stepping back" point is placed
        % second to last, so it's in the proper step order, like Fig. 3.7
        % in the notes.
        
        [M,~] = size(f_alpha);
        
        f_alpha_temp1 = f_alpha(end,:);
        f_alpha_temp2 = f_alpha(end-1,:);
        f_alpha(M,:) = f_alpha_temp2;
        f_alpha(M-1,:) = f_alpha_temp1;        
        
        alpha1 = f_alpha(1,2);
        alpha2 = f_alpha(2,2);
        alpha3 = f_alpha(3,2);
        
        f1 = f_alpha(1,1);
        f2 = f_alpha(2,1);
        f3 = f_alpha(3,1);
        
        delta_alpha = alpha2-alpha1;
        
        alpha_star = alpha2 + (delta_alpha*(f1 - f3))/(2*(f1 - 2*f2 + f3));
     end
 
    function [xnew,fnew] = take_step(obj,x,alpha_star,s)
        xnew = x + alpha_star*s;
        fnew = obj(xnew);
    end
    
    function [xnew,fnew] = Newton_quad(obj,gradobj,x)
        % REWRITE USING HESSIAN FUNCTION AND SUBS() TO SOLVE FOR THE VALUES
        % OF THE HESSIAN AT THE REQUESTED X POINTS
%         H = hessian(obj(x))
        syms x1 x2 x3
        % Enter the desired function (must match the function in
        % fminunDrivHW.m):
        % Quadratic test function
        f = 20 + 3*x1 - 6*x2 + 8*x3 + 6*x1^2 - 2*x1*x2 - x1*x3 + x2^2 + 0.5*x3^2;
        H = hessian(f,[x1,x2,x3]);
        H = double(H);
        
%         % Rosenbrock's function
%         f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
%         H = hessian(f,[x1,x2]);
%         H = double(H);
        
        grad = gradobj(x);
        delta_x = -inv(H)*grad;
        xnew = x + delta_x;
        fnew = obj(xnew);
    end