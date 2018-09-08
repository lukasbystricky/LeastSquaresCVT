function [generators1, mass] = cvt_step_lloyds_2d(generators0, f)

addpath('triangle_dunavant_rule');
tol = 1e-12;
[n,~] = size(generators0);
generators1 = zeros(size(generators0));

%% construct Voronoi diagram, find all polygons
v = voronoi(generators0(:,1),generators0(:,2));
[va,~] = voronoin(generators0);

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
    
    distances = sum((generators0 - vertice_rep).^2, 2);
    
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
    angles = atan2(vertices_poly(2,:) - generators0(g,2), ...
                vertices_poly(1,:) - generators0(g,1));
    
    [~, I] = sort(angles);
    
    for i = 1:length(angles)
        vertices_reordered(i) = vertices_ind(I(i));
    end
    
    polygons{g} = vertices_reordered;
end

%% integrate x*f(x) over each polygon, move generators 

mass = zeros(n,1);

for g = 1:n
    
    % split polygon into triangles
    region_vertices = vertices(:, polygons{g});
    vx = region_vertices(1,:);
    vy = region_vertices(2,:);
    
    cx = generators0(g,1);
    cy = generators0(g,2);
    
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
    
    center_of_mass = [0,0];
    mass_modified = 0;
    % integrate over each triangle
    
    rule = 20;
    
    for t = 1:n_vert
        v1 = [triangles(t,1,1),triangles(t,1,2)];
        v2 = [triangles(t,2,1),triangles(t,2,2)];
        v3 = [triangles(t,3,1),triangles(t,3,2)];
        
        [center_of_mass_tri, mass_tri, mass_modified_tri] = ...
                        integrate_triangle(v1, v2, v3, f, rule);
        
        center_of_mass = center_of_mass + center_of_mass_tri;
        mass_modified = mass_modified + mass_modified_tri;
        mass(g) = mass(g) + mass_tri;
        
    end
    
    generators1(g,:) = center_of_mass / mass_modified;
end



