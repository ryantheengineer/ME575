function [Lgrad] = Lagrangian_gradient(grad_f,lamda,grad_g)
    Lgrad = grad_f - lamda*grad_g;
end