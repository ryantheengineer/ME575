function [grad] = central_partial(fun, x)

    n = numel(x);
%     fb = fun(x);
    h = 0.0001;
    grad = zeros(n,1);
    for i = 1:n
        xtplus = x;
        xtminus = x;
        xtplus(i) = xtplus(i) + h;
        xtminus(i) = xtminus(i) - h;
        grad(i) = (fun(xtplus) - fun(xtminus))/(2*h);
    end
end