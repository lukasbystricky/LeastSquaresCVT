function [generators1, mass, energy] = cvt_step_lloyds_2d(generators0, rho)
%% cvt_step_lloyds_2d performs a single step of Lloyds algorithm
%
% Performs a single step of Lloyd's algorithm to update the Voronoi
% generators. This is accomplished using high order accurate quadrature
% rules.
%
% PARAMETERS:
% INPUTS:
% generators0, REAL (n,2), the initial Voronoi centroids
% rho, FUNCTION (x,y), the density function used to compute the mass
%
% OUTPUT:
% generators1, REAL (n,2), the updated Voronoi centroids
% mass, REAL (n,1), the mass of each region associated with generators0
% energy, REAL, the total energy associated with the centroids generators0

[n,~] = size(generators0);
generators1 = zeros(size(generators0));

[mass, center_of_mass, energy] = vornoi_compute_mass(generators0, rho);

for g = 1:n
    generators1(g,:) = center_of_mass(g, :);
end



