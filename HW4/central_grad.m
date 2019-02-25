function [constraint_partials] = central_grad(ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
    
    x = Elem(:,3);
    n = numel(x);
    h = 0.001;
    constraint_partials = zeros(n,n);
    
    % Repurpose forward_partial, which goes down the column of Elem, to do
    % this for each constraint
    
    % For each constraint:
    for j = 1:n
        % For each variable:
        for i = 1:n
            xtplus = x;
            xtminus = x;
            xtplus(i) = xtplus(i) + h;
            xtminus(i) = xtminus(i) - h;

            % Get the constraint values at the positive perturbation point
            for k=1:nelem
                Elem(k,3) = xtplus(k);
            end
            [~,cxtplus] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
            
            % Get the constraint values at the negative perturbation point
            for k=1:nelem
                Elem(k,3) = xtminus(k);
            end
            [~,cxtminus] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
            
            constraint_partials(i,j) = (cxtplus(j) - cxtminus(j))/(2*h);
        end
    end
    
end