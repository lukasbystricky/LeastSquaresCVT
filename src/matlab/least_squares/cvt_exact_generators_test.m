%% cvt_exact_generators_test
% Driver to create an exact CVT of the Christoffel function in 2D. 

addpath('..\exact_cvt\')
addpath('..\exact_cvt\triangle_dunavant_rule')

total_degree = 5;
n_generators = 50;
n_iterations = 100;
n_trials = 5;

condition_number_array = zeros(n_iterations, n_trials);
energy_array = zeros(n_iterations, n_trials);

M = compute_christoffel_m(total_degree);

Phi = struct();
Phi.index_set = total_degree_set(total_degree);
Phi.basis_card = size(Phi.index_set, 1);
Phi.value = @(X) eval_phi_fast(X, Phi.index_set);

for k = 1:n_trials
    generators0 = rejection_sampling_uniform(n_generators, Phi, false);
    
    [generators_exact, masses_exact, energy_exact] = cvt_lloyds_2d(generators0, ...
        n_iterations, @(x) christoffel(total_degree, x), true);
    
    energy_array(:, k) = energy_exact;
    
    condition_number_exact = zeros(n_iterations, 1);
    Phi = struct();
    Phi.index_set = total_degree_set(total_degree);
    Phi.basis_card = size(Phi.index_set, 1);
    Phi.value = @(X) eval_phi_fast(X, Phi.index_set);
    
    for i = 1:n_iterations
        G = create_gramian( generators_exact(:, :, i), Phi );
        condition_number_exact(i) = cond(G);
        condition_number_array(i, k) = cond(G);
    end
    
    save(['generators_exact_', num2str(n_generators), '_trial', num2str(k)],'generators_exact',...
        'masses_exact','energy_exact', 'condition_number_exact');
end