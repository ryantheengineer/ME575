function [f,gradf] = objfungrad_forward(Truss,forward_partial,Data)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    Data;

    % % insert areas (design variables) into correct matrix
    % for i=1:nelem
    %     Elem(i,3) = x(i);
    % end

    % call Truss to get weight and stresses
    [weight,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

    %objective function
    f = weight; %minimize weight

    % Gradient of the objective function
    if nargout > 1
        gradf = forward_partial(Truss,ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
    end
    
end

