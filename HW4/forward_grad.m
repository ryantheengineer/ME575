function [constraint_partials] = forward_grad(ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
    
    x = Elem(:,3);
    n = numel(x);
    [~,cb] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
    h = 0.01;
    constraint_partials = zeros(n,n);
    
    % Repurpose forward_partial, which goes down the column of Elem, to do
    % this for each constraint
    
    % For each constraint:
    for j = 1:n
        % For each variable:
        for i = 1:n
            xt = x;
            xt(i) = xt(i) + h;

            for k=1:nelem
                Elem(k,3) = xt(k);
            end
            [~,cxt] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
            constraint_partials(i,j) = (cxt(j) - cb(j))/h;
        end
    end
    
end