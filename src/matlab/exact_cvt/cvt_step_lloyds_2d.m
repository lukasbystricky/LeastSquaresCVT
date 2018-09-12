function [generators1, mass, energy] = cvt_step_lloyds_2d(generators0, f)

[n,~] = size(generators0);
generators1 = zeros(size(generators0));

[mass, center_of_mass, energy] = vornoi_compute_mass(generators0, f);

for g = 1:n
    generators1(g,:) = center_of_mass(g, :);
end



