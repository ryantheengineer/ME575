function [update] = BFGS_Hessian(gamma_k,H_k,delta_xk)
term2 = (gamma_k*transpose(gamma_k))/(transpose(gamma_k)*delta_xk);

term3 = (H_k*delta_xk*transpose(delta_xk)*H_k)/(transpose(delta_xk)*H_k*delta_xk);

update = H_k + term2 - term3;


end