%% cvt_kmeans_generator_test
% Driver to test the accuracy of the kmeans sampling algorithm for the 
% Christoffel function in 2D. The energy and the condition number of the 
% simulations are compared to an exact CVT. 

addpath('..\exact_cvt\')
addpath('..\exact_cvt\triangle_dunavant_rule')

total_degree = 5;
n_iterations = 100;
sample_type = 'cvt';

load('generators_exact_50_trial1');
generators0 = generators_exact(:,:,1);
[n_generators, ~] = size(generators0);

Phi = struct();
Phi.index_set = total_degree_set(total_degree);
Phi.basis_card = size(Phi.index_set, 1);
Phi.value = @(X) eval_phi_fast(X, Phi.index_set);
n_samples = 100*n_generators;

n_trials = 50;

energy_kmeans_array = zeros(n_iterations, n_trials);
condition_number_kmeans_array = zeros(n_iterations, n_trials);

for k = 1 : n_trials
    [generators_kmeans, energy_kmeans_array(:,k), masses_kmeans] = cvt_kmeans(...
                    generators0, Phi, n_samples, sample_type, n_iterations, false);

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
xlabel('iteration');
title('Energy')

matlab2tikz(['energy_generators_', num2str(n_generators), '_', ...
                    sample_type, '_samples_', num2str(n_samples/n_generators),...
                    '.tex'], 'standalone', true, 'floatFormat', '%.6f');
                
figure()
plot(condition_number_kmeans_array, 'b');
hold on
plot(mean(condition_number_kmeans_array, 2), 'b', 'Linewidth', 2);
plot(condition_number_exact, 'r', 'Linewidth', 2);
title('cond(G)')
xlabel('iteration');

matlab2tikz(['condition_number_generators_', num2str(n_generators), '_', ...
                    sample_type, '_samples_', num2str(n_samples/n_generators),...
                    '.tex'], 'standalone', true, 'floatFormat', '%.6f');