clear all;
close all;

% Choose a starting design
xs = [5;5];
fs = objective(xs);

% Select Ps, Pf, N, and calculate Ts, Tf, and F
Ps = 0.8;               % Probability of acceptance at start
Pf = 0.001;             % Probability of acceptance at finish
N = 100;                % Number of cycles

Ts = -1/log(Ps);        % Temperature at start
Tf = -1/log(Pf);        % Temperature at finish
F = (Tf/Ts)^(1/(N-1));  % Temperature reduction factor each cycle

% Perturbation information
delta = 0.5;           % Max perturbation
n = 5;                  % Starting number of perturbations per cycle

% Holding variables
dE = 0;
dElast = 0;
dEavg = 0;
cycles = [1:N];
objvals = zeros(1,N);
Pvec = [0];

% Set starting values
xc = xs;
fc = fs;
T = Ts;

% Step through the cycles
for i = 1:N
    % Add the current objective value to the objective vector for plotting
    objvals(i) = objective(xc);
    
    % Step through the perturbations
    for j = 1:n
        % Perturb xc by some random value within delta
        x1p = unifrnd(-delta,delta);
        x2p = unifrnd(-delta,delta);
        perturb = [x1p; x2p];
        xp = xc + perturb;
        
        % Get the objective value at the perturbed point
        fp = objective(xp);
        
        % Calculate values for Boltzmann function in case they're needed
        dE = abs(fp-fc);
        dElast = dE;
        if i == 1 && j == 1
            dEavg = dE;
        else
            dEavg = (dEavg + dElast)/2;
        end
        
        P = Boltzmann(dE,dEavg,T);
        % Save the value of P
        Pvec = [Pvec;P];
        
        % Check if the new design is better than the old design
        if fp < fc
            xc = xp;    % Accept as current design if better
        else
            % If the new design is worse, generate a random number and
            % compare to the Boltzmann probability. If the random number is
            % lower than the Boltzmann probability, accept the worse design
            % as the current design
            randnum = unifrnd(0,1);
            
            if randnum < P
                xc = xp;
            end
        end     
    end
    % Decrease the temperature by factor F
    T = F*T;
    % Increase the number of perturbations each cycle
    n = n+2;
    % Decrease the maximum perturbation each cycle
    delta = delta*0.9;
end
        

% Plot the cooling curve
figure(1);
plot(cycles,objvals);
xlabel('Cycles');
ylabel('Objective');

% figure(2);
% plot(Pvec(2:end));






% % Check the function by creating a 3D plot
% [X1,X2] = meshgrid(-5:0.1:5,-5:0.1:5);
% F = 2 + 0.2*X1.^2 + 0.2*X2.^2 - cos(pi*X1) - cos(pi*X2);
% surf(X1,X2,F)




%% Functions

function [f] = objective(x)
    f = 2 + 0.2*x(1).^2 + 0.2*x(2).^2 - cos(pi.*x(1)) - cos(pi.*x(2));
end

% Boltzmann probability
function [P] = Boltzmann(dE,dEavg,T)
    P = exp(-dE/(dEavg*T));
end