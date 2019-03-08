%% HW5.m:
% Simulated annealing algorithm

% Algorithm steps, from 5.5.2.3:
% 1. Choose a starting design
xc = [5;5];  % ARBITRARY RIGHT NOW but must be within (-5,5)
fc = objective(xc);
n = 10;     % number of perturbations

% 2. Select Ps, Pf, N, and calculate Ts, Tf, and F
Ps = 0.7;
Pf = 0.001;
N = 100;     % number of cycles

delta = 0.5;    % Is this the max perturbation?

Ts = -1/log(Ps);
Tf = -1/log(Pf);
F = (Tf/Ts)^(1/(N-1));

% 3. Randomly perturb the design to different discrete values close to the
% current design
dElast = 0;

cycles = [1:N];
objvals = zeros(1,N);

for i = 1:N
    % Add the current objective value to the objective vector
    objvals(i) = objective(xc);
    for j = 1:n
        % Randomly perturb x within limits set by delta
        perturb = [unifrnd(-delta,delta); unifrnd(-delta,delta)];
        xp = xc + perturb;
        
        % Check if new design is better than the old design
        fp = objective(xp);
        if fp < fc
            % 4. If the new design is better, accept it as the current design
            xc = xp;
        else
            dE = fp - fc;   % MIGHT NEED TO CHANGE THE ORDER, CHECK LATER
            if i == 1 && j == 1
                dEavg = dE;
                T = Ts;
            else
                dEavg = (dEavg + dElast)/2;
            end
            
            % 5. If the new design is worse, generate a random number 
            % between 0 and 1 using a uniform distribution. Compare this 
            % number to the Boltzmann probability. If the random number is
            % lower than the Boltzmann probability, accept the worse design
            % as the current design.
            randnum = unifrnd(0,1);
            P = Boltzmann(dE,dEavg,T);     % Get Boltzmann probability
            if randnum < P
                xc = xp;                
            end
        end
       dElast = dE; 
    end
    % 7. Decrease temperature according to T(n+1) = F*Tn
    T = F*T;
end


% Plot the cooling curve
plot(cycles,objvals);
% As a check, plot the probabilities, F, and other values that should
% decrease or have some trend over time here (I think I forgot to update F
% after each cycle)


%% Functions
% Objective function
function [f] = objective(x)
    f = 2 + 0.2*x(1)^2 + 0.2*x(2)^2 - cos(pi*x(1)) - cos(pi*x(2));
end

% Boltzmann probability
function [P] = Boltzmann(dE,dEavg,T)
    P = exp(-dE/(dEavg*T));
end