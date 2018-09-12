addpath('..\exact_cvt\')
addpath('..\exact_cvt\triangle_dunavant_rule')

total_degree = 5;
n_generators = 50;
n_iterations = 100;

M = compute_christoffel_m(total_degree);
generators0 = initial_generators(n_generators, M, ...
                @(x,y) christoffel(total_degree, x, y));
 
[generators_exact50, masses_exact50] = cvt_lloyds_2d(generators0, ...
            n_iterations, @(x,y) christoffel(total_degree, x, y));
        
        
condition_number = zeros(n_iterations, 1);
Phi = struct();
Phi.index_set = total_degree_set(total_degree);
Phi.basis_card = size(Phi.index_set, 1);
Phi.value = @(X) eval_phi_fast(X, Phi.index_set);

for i = 1:n_iterations
    G = create_gramian( generators_exact50(:, :, i), Phi );
    condition_number(i) = cond(G);
end