function [c,ceq,DC,DCeq] = confungrad_forward(Truss,forward_partial,Data)

    Data;
    
    % call Truss to get stresses (constraints)
    [~,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
    
    c = zeros(nelem,1);
    for i = 1:nelem
        c(i) = stress(i);
    end
    
    % No nonlinear equality constraints
    ceq = [];

    % Gradient of the constraints
    if nargout > 2
        DC = forward_grad(Truss,ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
        DCeq = [];
    end
    
end