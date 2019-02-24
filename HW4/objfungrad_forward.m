function [f,gradf] = objfungrad_forward(Truss,forward_partial,Data)

    Data;

    % % insert areas (design variables) into correct matrix
    % for i=1:nelem
    %     Elem(i,3) = x(i);
    % end

    % call Truss to get weight (objective)
    [weight,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

    %objective function
    f = weight; %minimize weight

    % Gradient of the objective function
    if nargout > 1
        gradf = forward_partial(Truss,ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
    end
    
end

