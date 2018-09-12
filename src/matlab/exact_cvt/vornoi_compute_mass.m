function [mass, center_of_mass, energy] = vornoi_compute_mass(generators, rho)
%% voronoi_compute_mass computes the masses and centers of each Voronoi region
%
% Computes the masses and centers of masses of each Voronoi region, as well
% as the total energy associated with the generators. The mass of region I
% is \int_I \rho(x,y) dx dy. The center of mass  is 
% \int_I x_i \rho(x,y) dx dy / \int_I \rho(x,y) dx dy. The total energy is
% \sum_I (\int_I \rho(x,y) || (x,y) - *cx_I, cy_I||^2 dx dy). These integrals
% are computed using the Dunavant quadrature rules that can integrate all
% polynomials up to total degree 20 exactly.
%
% PARAMETERS:
% INPUTS:
% generators, REAL (n,2), the Voronoi generators, i.e. the centroids of
%   each region
% rho, FUNCTION (x,y), the density function used to compute the mass
%
% OUTPUT:
% mass, REAL (n), the masses of each Voronoi region
% center_of_mass, REAL (n,2), the center of mass of each Voronoi region
% energy, REAL, the total energy

tol = 1e-12;
[n,~] = size(generators);

%% construct Voronoi diagram, find all polygons
v = voronoi(generators(:,1),generators(:,2));
[va,~] = voronoin(generators);

% find all vertices of Voronoi regions, including on unit box
vL = [v(2).XData; v(2).YData];
bL = [-1, 1, 1, -1, -1; -1, -1, 1, 1, -1];

% find intersection of Voronoi boundaries with unit box
intersect = InterX(vL, bL);
vertices = [intersect, [-1, -1, 1, 1; -1, 1, 1, -1]];

% remove vertices outside unit box
va = va(va(:,1) > -1 & va(:,1) < 1 & va(:,2) > -1 & va(:,2) < 1,:);
vertices = [vertices, va'];

% find distance from generators to each vertice
polygons = cell(n,1);
for i = 1:size(vertices,2)
    vertice_rep = repmat(vertices(:,i)', n, 1);
    
    distances = sum((generators - vertice_rep).^2, 2);
    
    distances = distances - min(distances);
    
    for g = 1:n
        if distances(g) < tol
            polygons{g} = [polygons{g},i];
        end
    end
end

% reorder polygon vertices to get convex polygon
for g = 1:n
    vertices_poly = vertices(:, polygons{g});
    vertices_ind = polygons{g};
    vertices_reordered = zeros(size(polygons{g}));
    angles = atan2(vertices_poly(2,:) - generators(g,2), ...
                vertices_poly(1,:) - generators(g,1));
    
    [~, I] = sort(angles);
    
    for i = 1:length(angles)
        vertices_reordered(i) = vertices_ind(I(i));
    end
    
    polygons{g} = vertices_reordered;
end

%% integrate x*f(x) over each polygon, move generators 

mass = zeros(n,1);
center_of_mass = zeros(n,2);

energy = 0;
for g = 1:n
    
    % split polygon into triangles
    region_vertices = vertices(:, polygons{g});
    vx = region_vertices(1,:);
    vy = region_vertices(2,:);
    
    cx = generators(g,1);
    cy = generators(g,2);
    
    n_vert = length(vx);
    triangles = zeros(n,3,2);
    
    for i = 1:n_vert
        triangles(i,1,1) = cx;
        triangles(i,1,2) = cy;
        triangles(i,2,1) = vx(i);
        triangles(i,2,2) = vy(i);
        
        if i + 1 > n_vert
            triangles(i,3,1) = vx(1);
            triangles(i,3,2) = vy(1);
        else
            triangles(i,3,1) = vx(i+1);
            triangles(i,3,2) = vy(i+1);
        end
    end
    
    mass_modified = 0;
    % integrate over each triangle
    
    rule = 20;
    
    for t = 1:n_vert
        v1 = [triangles(t,1,1),triangles(t,1,2)];
        v2 = [triangles(t,2,1),triangles(t,2,2)];
        v3 = [triangles(t,3,1),triangles(t,3,2)];
        
        [center_of_mass_tri, mass_tri, mass_modified_tri, energy_tri] = ...
                        integrate_triangle(v1, v2, v3, rho, rule);
        
        center_of_mass(g,:) = center_of_mass(g,:) + center_of_mass_tri;
        mass_modified = mass_modified + mass_modified_tri;
        mass(g) = mass(g) + mass_tri;
        energy = energy + energy_tri;
    end
    
    center_of_mass(g,:) = center_of_mass(g,:) / mass_modified;
end
