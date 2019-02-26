function [grad] = forward_diff(x)
% Testing the fundamental algorithm behind forward difference method as a
% gut check

    n = numel(x);
    fxb = sin(x);
    h = 0.001;
    grad = zeros(n,1);
    for i = 1:n
        % temporary holding variable for perturbation
        xt = x;
        xt(i) = xt(i) + h;
        fxt = sin(xt);
        grad(i) = (fxt - fxb)/h;
    end
    
end