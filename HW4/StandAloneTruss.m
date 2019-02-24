% program to find displacements and stresses in a truss us FEM
clear;
Data;

% x = Elem(:,3);

[weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
stress
weight

% Testing derivatives
grad = forward_partial(@Truss,ndof, nbc, nelem, E, dens, Node, force, bc, Elem)