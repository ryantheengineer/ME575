clear;
x1 = [-1.695; 2.157];

x2 = [-0.592; -1.162];
f2 = objective(x2)
g2 = constraint(x2)
s2 = 0.230
lamda = 0

gradf1 = gradient_f(x1)
gradg1 = gradient_g(x1)

gradf2 = gradient_f(x2)
gradg2 = gradient_g(x2)



Lgrad1 = Lagrangian_gradient(gradf1,lamda,gradg1)
Lgrad2 = Lagrangian_gradient(gradf2,lamda,gradg2)

gamma1 = Lgrad2 - Lgrad1

deltax1 = x2 - x1

LagrHess1 = [20.762, 5.629;
             5.629, 1.910];

LagrHess2 = BFGS_Hessian(gamma1,LagrHess1,deltax1)

Jacobian = [transpose(gradg2)]

% Assemble the coefficient matrix:
coeffMat = [LagrHess2(1,1),LagrHess2(1,2),0,-Jacobian(1);
            LagrHess2(2,1),LagrHess2(2,2),0,-Jacobian(2);
            0, 0, lamda, s2;
            Jacobian(1),Jacobian(2), -1, 0]
        
        
% Calculate the residual vector
mu = 5/25;
resid2 = zeros(4,1);
resid2(1) = gradf2(1) - lamda*gradg2(1);
resid2(2) = gradf2(2) - lamda*gradg2(2);
resid2(3) = s2*lamda - mu;
resid2(4) = g2 - s2;
resid2

% Solve for the delta vector using coeffMat and resid2
delta_vec = coeffMat\(-resid2)

x3 = x2 + delta_vec(1:2)
s3 = s2 + delta_vec(3)
lamda3 = lamda + delta_vec(4)

f3 = objective(x3)
g3 = constraint(x3)

P = f3 + 2*abs(g3)


%% Functions
function [f] = objective(x)
    f = x(1)^4 - 2*x(2)*x(1)^2 + x(2)^2 + x(1)^2 - 2*x(1) + 5;
end

function [g] = constraint(x)
    g = -(x(1) + 0.25)^2 + 0.75*x(2);
end

function [gradf] = gradient_f(x)
    gradf = [4*x(1)^3 - 4*x(1)*x(2) + 2*x(1) - 2;
             -2*x(1)^2 + 2*x(2)];
end

function [gradg] = gradient_g(x)
    gradg = [(-2*x(1)-0.5);
             0.75];
end

function [Lgrad] = Lagrangian_gradient(grad_f,lamda,grad_g)
    Lgrad = grad_f - lamda*grad_g;
end

function [update] = BFGS_Hessian(gamma_k,H_k,delta_xk)
    term2 = (gamma_k*transpose(gamma_k))/(transpose(gamma_k)*delta_xk);

    term3 = (H_k*delta_xk*transpose(delta_xk)*H_k)/(transpose(delta_xk)*H_k*delta_xk);

    update = H_k + term2 - term3;
end