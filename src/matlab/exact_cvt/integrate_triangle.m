function [center_of_mass, mass, modified_mass, energy] = ...
                    integrate_triangle(v1, v2, v3, rho, rule)
%% integrate_triangle computes the mass and center of mass of a triangle
%
% Computes the mass, center of mass, Voronoi energy and the CVT modified
% mass of a trianglular part of a Voronoi cell. The mass of triangle T
% is \int_T \rho(x,y) dx dy and its center of mass  is 
% \int_T x_i \rho(x,y) dx dy / \int_T \rho(x,y) dx dy. The modified mass of
% T is \int_T \rho(x,y)^2 dx dy. The energy is
% \int_T \rho(x,y) || (x,y) - (v1_x, v1_y)||^2 dx dy. These integrals
% are computed using the Dunavant quadrature rules that can integrate all
% polynomials up to total degree 20 exactly.
%
% PARAMETERS:
% INPUTS:
% v1, REAL (2), a vertice of the triangular region, v1 is also the center
%   of the Voronoi region
% v2, REAL (2), the second vertice of the triangular region
% v3, REAL (2), the third vertice of the triangular region
% rho, FUNCTION (x,y), the density function which determines the mass 
% rule, INTEGER, the order accuracy of the quadrature rule
%
% OUTPUT:

% center_of_mass, REAL (2), the center of mass of the triangle
% mass, REAL, the mass of the triangle
% modified_mass, REAL, the modified mass of the triangle
% energy, REAL, the energy of the triangle

order_num = dunavant_order_num ( rule );
[ xy, w ] = dunavant_rule ( rule, order_num );

% map to physical triangle
xy2 = reference_to_physical_t3 ( [v1;v2;v3]', order_num, xy );
area = triangle_area ( [v1; v2; v3]' );

% perform quadrature
f_quad = rho([xy2(1,:); xy2(2,:)]')';

center_of_mass(1) = sum(w .* f_quad.^2 .* xy2(1,:));
center_of_mass(2) = sum(w .* f_quad.^2 .* xy2(2,:));
modified_mass = sum(w .* f_quad.^2);
mass = sum(w .* f_quad);
energy = sum(w .* f_quad.^2 .* sum( (xy2' - repmat( v1, size(xy2, 2), 1)).^2, 2)');
%v1 is always the center of the Voronoi region

mass = area * mass;
modified_mass = area * modified_mass;
center_of_mass = area * center_of_mass;
energy = area * energy;