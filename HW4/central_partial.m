function [grad] = central_partial(Truss,ndof, nbc, nelem, E, dens, Node, force, bc, Elem)

    x = Elem(:,3);
    n = numel(x);
    h = 0.0001;
    grad = zeros(n,1);
    for i = 1:n
        xtplus = x;
        xtminus = x;
        xtplus(i) = xtplus(i) + h;
%         xtplus
        xtminus(i) = xtminus(i) - h;
%         xtminus
        
        % Get the objective value at the positive perturbation point
        for j=1:nelem
            Elem(j,3) = xtplus(j);
        end
        [fxtplus,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
        
        % Get the objective value at the negative perturbation point
        for j=1:nelem
            Elem(j,3) = xtminus(j);
        end
        [fxtminus,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
        
        grad(i) = (fxtplus - fxtminus)/(2*h);
    end
end