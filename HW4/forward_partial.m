function [grad] = forward_partial(Truss,ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
    
    x = Elem(:,3);
    n = numel(x);
    [fb,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
    h = 0.0001;
    grad = zeros(n,1);
    for i = 1:n
        xt = x;
        xt(i) = xt(i) + h;
%         xt
        for j=1:nelem
            Elem(j,3) = xt(j);
        end
        [fxt,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
        grad(i) = (fxt - fb)/h;
    end
    
end