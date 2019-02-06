
     function [xopt, fopt, exitflag] = fminun(obj, gradobj, x0, stoptol, algoflag)
        global nobj
        % get function and gradient at starting point
        [n,~] = size(x0); % get number of variables
        f = obj(x0)
        grad = gradobj(x0)
        x = x0;
        searchcount = 0;
        
        %set starting step length
        alpha = 0.05;
     
        if (algoflag == 1)      % steepest descent
            while nobj < 20000
%                 disp('New search direction');
%                 searchcount = searchcount + 1
                s = srchsd(grad)
                alpha_star = minimizing_step(obj,s,x,f,alpha);
                [xnew,fnew] = take_step(obj,x,alpha_star,s);
                % Check to see if the gradient is within the desired tolerance
                grad = gradobj(xnew)
                if abs(grad(1))>stoptol || abs(grad(2))>stoptol || abs(grad(3))>stoptol % if the gradient is not smaller than stoptol
                    x = xnew;
                    f = fnew;
                    alpha = 0.05;
                    if searchcount > 2000
                        exitflag = 0;
                        break
                    end
                    continue
                else
                    exitflag = 1;
                    break
                end
            end
            
        end
        
%         if (algoflag == 2)      % BFGS quasi-Newton
%             s = srchbfgs(grad)
%         end
        
%             if all(abs(grad)>stoptol)   % if the gradient is not smaller than the

        % while the number of function calls is less than our set maximum
%         while nobj < 1000   
            
            
            
                
        
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
                count = count + 1;
                continue
            else
                % if the function begins to increase, then curve fit the
                % data with a parabola, and step to the minimum of the
                % parabola.
                alpha = alpha/2;
                xnew = x +alpha*s
                fnew = obj(xnew);
                f_alpha_temp = [fnew,alpha];
                f_alpha = [f_alpha; f_alpha_temp];
                
%                 % (CHECK IF THE LAST POINT IS IN BETWEEN THE PREVIOUS TWO
%                 % POINTS)
%                 if x_f_alpha(end,2)
%                 
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
        
        % Once we have a set of points that bracket the true minimum, we
        % pick the minimum of our points and the two that are next to it in
        % the line search (how to do this I'm not sure yet). Then we use
        % the formula found on page 18 of chapter 3 in the notes to solve
        % for alpha_star
        
        
        
        
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