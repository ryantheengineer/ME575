% program to find displacements and stresses in a truss us FEM
clear;
Data;

% x = Elem(:,3);

% Elem(1,3) = 10;

[weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
stress
weight

% Testing derivative of objective
% grad_forward = forward_partial(@Truss, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
% grad_central = central_partial(@Truss, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)

% Testing gradient of constraints
constraints_forward = forward_grad(@Truss, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
constraints_central = central_grad(@Truss, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)