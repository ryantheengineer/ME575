
     function [xopt, fopt, exitflag] = fminun(obj, gradobj, x0, stoptol, algoflag)
        global nobj
        global ngrad
        % get function and gradient at starting point
        [n,~] = size(x0); % get number of variables
        f = obj(x0)
        grad = gradobj(x0)
        x = x0;
        xsearch = [x0];
     
        if (algoflag == 1)      % steepest descent
            while ngrad < 100
                % Set starting step length
%                 alpha = 0.12; % step length for quadratic
                alpha = 0.01 % step length for Rosenbrock's

                s = srchsd(grad);
                alpha_star = linesearch(obj,s,x,f,alpha);
                [xnew,fnew] = take_step(obj,x,alpha_star,s);
                xsearch = [xsearch,xnew];   % save x vectors that define search pattern
                gradnew = gradobj(xnew);
                
%                 %Check if alpha_star is correct (should get zero here)
%                 alpha_star_check = s'*gradnew;
                
                % Check to see if the gradient is within the desired tolerance
                grad = gradobj(xnew)
                x = xnew;
                f = fnew;
                
                if all(abs(grad) > stoptol)
                    continue
                else
                    exitflag = 1;
                    break
                end
            end
        end
        
        % BFGS quasi-Newton
        if (algoflag == 2)
            % Define initial N and search direction (s normalized)
            N = eye(n);
            s = srchbfgs(grad,N);
            
            while ngrad < 2000
                % Execute line search in direction s, and get xnew and gradnew
%                 alpha = 0.12;      % alpha for quadratic function
                alpha = 0.001;     % alpha for Rosenbrock's function
                alpha_star = linesearch(obj,s,x,f,alpha);
                [xnew,fnew] = take_step(obj,x,alpha_star,s);
                xsearch = [xsearch,xnew];   % save x vectors that define search pattern
                gradnew = gradobj(xnew);
                
                % Check if alpha_star is correct (should get zero here)
                alpha_star_check = s'*gradnew;
                
                % Solve for delta_x and gamma
                delta_x = xdist(x,xnew);
                gamma = get_gamma(grad,gradnew);

                % Solve for Nnew
                Nnew = BFGS_update(N,gamma,delta_x);

                % Solve for new s and repeat
                grad = gradnew
                x = xnew;
                f = fnew;
                N = Nnew;
                s = srchbfgs(grad,N);
                
                % Check to see if gradnew is within the desired tolerance
                if all(abs(grad) > stoptol)
                    continue
                else
                    exitflag = 1;
                    break
                end
            end
        end            

        save Rosenbrock_search.mat xsearch
        
        xopt = xnew;
        fopt = fnew;
        
        exitflag = 1;   % This is automatically set to stop the solver from iterating for a long time, and it's why xopt doesn't result in a minimum right now
     end
     
     
     %% Steepest descent functions
     % get steepest descent search direction as a column vector
     function [s] = srchsd(grad)
        mag = sqrt(grad'*grad);
        s = -grad/mag;
     end
     
     function [alpha_star] = linesearch(obj,s,x,f,alpha)
        % Create holding vectors for the f and alpha values
        alpha_vector = zeros(3,1);
        f_vector = zeros(3,1);
        
        alpha_prev = alpha;
        f_prev = f;
        
        count = 0;
        
        while count < 1000   % this might not be the right limit
            % take a step
            [~,fnew] = take_step(obj,x,alpha,s);
            
            % check if the value of fnew is less than fprevstep (is the function
            % still decreasing?):
            if fnew < f_prev
                % Save the latest values for fnew as the first end of the
                % bracket around the minimum
                alpha_vector(1) = alpha;
                f_vector(1) = fnew;
                
                alpha_prev = alpha;
                f_prev = f;
                
                % Set the next alpha
                alpha = 2*alpha;
                
                count = count + 1;
                
                continue
            else
                % Set the far end of the bracketing vectors
                alpha_vector(3) = alpha;
                f_vector(3) = fnew;
                
                % Take a new step that is the average of the previous two
                % steps
                alpha = (alpha + alpha_prev)/2;
                [~,fnew] = take_step(obj,x,alpha,s);
                
                alpha_vector(2) = alpha;
                f_vector(2) = fnew;
                  
                break
            end
        end
        
        % Renaming for brevity
        alpha1 = alpha_vector(1);
        alpha2 = alpha_vector(2);
        alpha3 = alpha_vector(3);
        
        f1 = f_vector(1);
        f2 = f_vector(2);
        f3 = f_vector(3);
        
        delta_alpha = alpha2-alpha1;
        
        alpha_star = alpha2 + (delta_alpha*(f1 - f3))/(2*(f1 - 2*f2 + f3));
     end
 
    function [xnew,fnew] = take_step(obj,x,alpha_star,s)
        xnew = x + alpha_star*s;
        fnew = obj(xnew);
    end
        
    %% BFGS-specific functions
    % Get search direction using quasi-Newton method
    function [s] = srchbfgs(grad,N)
        s = -N*grad;
        mag = sqrt(s'*s);
        s = s/mag;
    end
    
    % Solve for gamma
    function [gamma] = get_gamma(grad,gradnew)
        gamma = gradnew-grad;
    end
    
    function [delta_x] = xdist(x,xnew)
        delta_x = xnew - x;
    end
    
    % BFGS update
    function [Nnew] = BFGS_update(N,gamma,delta_x)
        firstterm = 1 + (gamma'*N*gamma)/(delta_x'*gamma);
        secondterm = (delta_x*delta_x')/(delta_x'*gamma);
        thirdterm = (delta_x*gamma'*N + N*gamma*delta_x')/(delta_x'*gamma);
        
        Nnew = N + firstterm*secondterm - thirdterm;
    end