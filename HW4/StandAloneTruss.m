% program to find displacements and stresses in a truss us FEM
clear;
Data;

% for i=1:nelem
%     Elem(i,3) = 50;
% end

% Elem(1,3) = 7.94;
% Elem(2,3) = 0.1;
% Elem(3,3) = 8.06;
% Elem(4,3) = 3.94;
% Elem(5,3) = 0.1;
% Elem(6,3) = 0.1;
% Elem(7,3) = 5.74;
% Elem(8,3) = 5.57;
% Elem(9,3) = 5.57;
% Elem(10,3) = 0.1;

[weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
stress
weight

% Testing derivative of objective
grad_forward = forward_partial(ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
grad_central = central_partial(ndof, nbc, nelem, E, dens, Node, force, bc, Elem)

% Testing gradient of constraints
% constraints_forward = forward_grad(ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
% constraints_central = central_grad(ndof, nbc, nelem, E, dens, Node, force, bc, Elem)