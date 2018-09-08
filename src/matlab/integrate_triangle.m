function [center_of_mass, mass, modified_mass] = integrate_triangle(v1, v2, v3, f, rule)

% generate rule on reference triangle
order_num = dunavant_order_num ( rule );
[ xy, w ] = dunavant_rule ( rule, order_num );

% map to physical triangle
xy2 = reference_to_physical_t3 ( [v1;v2;v3]', order_num, xy );
area = triangle_area ( [v1; v2; v3]' );

% perform quadrature
center_of_mass = [0,0];
modified_mass = 0;
mass = 0;

for order = 1 : order_num
    
    x = xy2(1,order);
    y = xy2(2,order);
    
    center_of_mass(1) = center_of_mass(1) + w(order) * f(x,y).^2 * x;
    center_of_mass(2) = center_of_mass(2) + w(order) * f(x,y).^2 * y;
    modified_mass = modified_mass + w(order) * f(x,y).^2;
    mass = mass + w(order) * f(x,y);
    
end

mass = area * mass;
modified_mass = area * modified_mass;
center_of_mass = area * center_of_mass;
