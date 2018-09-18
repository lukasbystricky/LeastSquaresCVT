%% cvt_kmeans_generator_test
% Driver to test the accuracy of the kmeans sampling algorithm for the 
% Christoffel function in 2D. The energy and the condition number of the 
% simulations are compared to an exact CVT. 

addpath('..\exact_cvt\')
addpath('..\exact_cvt\triangle_dunavant_rule')

total_degree = 5;
n_iterations = 100;

load('generators_exact_50_trial1');
generators0 = generators_exact(:,:,1);
[n_generators, ~] = size(generators0);

Phi = struct();
Phi.index_set = total_degree_set(total_degree);
Phi.basis_card = size(Phi.index_set, 1);
Phi.value = @(X) eval_phi_fast(X, Phi.index_set);
num_sample = 1000*n_generators;

n_trials = 5;

energy_kmeans_array = zeros(n_iterations, n_trials);
condition_number_kmeans_array = zeros(n_iterations, n_trials);

for k = 1 : n_trials
    [generators_kmeans, energy_kmeans_array(:,k), masses_kmeans] = cvt_kmeans(...
                    generators0, Phi, num_sample, 'uniform', n_iterations, false);

    for i = 1:n_iterations
        G = create_gramian( generators_kmeans(:, :, i), Phi );
        condition_number_kmeans_array(i, k) = cond(G);
    end
    
end

figure()
plot(energy_kmeans_array, 'b');
hold on
plot(mean(energy_kmeans_array, 2), 'b', 'Linewidth', 2);
plot(energy_exact, 'r', 'Linewidth', 2);
title('Energy')

figure()
plot(condition_number_kmeans_array, 'b');
hold on
plot(mean(condition_number_kmeans_array, 2), 'b', 'Linewidth', 2);
plot(condition_number_exact, 'r', 'Linewidth', 2);
title('Condition number')