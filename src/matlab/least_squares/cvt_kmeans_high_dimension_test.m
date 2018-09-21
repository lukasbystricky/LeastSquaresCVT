%% cvt_kmeans_generator_test
% Driver to test the accuracy of the kmeans sampling algorithm for the 
% Christoffel function in 2D. The energy and the condition number of the 
% simulations are compared to an exact CVT. 

addpath('..\exact_cvt\')

total_degree = 5;
dimension = 2;
n_iterations = 100;
n_generators = 50;
sample_type = 'uniform';

Phi = struct();
Phi.index_set = total_degree_set(total_degree, dimension);
Phi.basis_card = size(Phi.index_set, 1);
Phi.value = @(X) eval_phi_fast(X, Phi.index_set);
n_samples = 100*n_generators;

n_trials = 1;

energy_kmeans_array = zeros(n_iterations, n_trials);
condition_number_kmeans_array = zeros(n_iterations, n_trials);


for k = 1 : n_trials
    
    
    generators0 = rejection_sampling_uniform(n_generators, Phi, false);
    [generators_kmeans, energy_kmeans_array(:,k), masses_kmeans] = cvt_kmeans(...
                    generators0, Phi, n_samples, sample_type, n_iterations,...
                    true, 'd2_td_5_g_50');

    for i = 1:n_iterations
        G = create_gramian( generators_kmeans(:, :, i), Phi );
        condition_number_kmeans_array(i, k) = cond(G);
    end
    
end