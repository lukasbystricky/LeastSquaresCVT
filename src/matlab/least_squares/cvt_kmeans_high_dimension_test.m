%% cvt_kmeans_generator_test
% Driver to test the accuracy of the kmeans sampling algorithm for the 
% Christoffel function in 2D. The energy and the condition number of the 
% simulations are compared to an exact CVT. 

addpath('..\exact_cvt\')

total_degree = 2;
dimension = 10;
n_iterations = 30;
n_generators = 100;
sample_type = 'uniform';

Phi = struct();
Phi.index_set = total_degree_set(total_degree, dimension);
Phi.basis_card = size(Phi.index_set, 1);
Phi.value = @(X) eval_phi_fast(X, Phi.index_set);
n_samples = [10, 100, 1000]*n_generators;

n_trials = 3;

energy_kmeans_array = zeros(n_iterations, n_trials);
condition_number_kmeans_array = zeros(n_iterations, n_trials);
masses_kmeans = zeros(n_generators, n_iterations, n_trials);

movements = zeros(n_iterations - 1, n_trials);

rng(1000); % set seed

generators0 = rejection_sampling_uniform(n_generators, Phi, false);

for k = 1 : n_trials
    
    [generators_kmeans, energy_kmeans_array(:,k), masses_kmeans(:, :,  k), areas_kmeans] = cvt_kmeans(...
                    generators0, Phi, n_samples(k), sample_type, n_iterations, ...
                    false, 'd3_td_5_g_200');

    for i = 1:n_iterations
        
        if i > 1
            d = generators_kmeans(:, :, i) - generators_kmeans(:, :, i - 1);
            movements(i - 1, k) = sqrt(sum(sum(d.^2, 2)));
        end
        
        
        G = create_gramian( generators_kmeans(:, :, i), Phi, [] );        
        condition_number_kmeans_array(i, k) = cond(G);
    end    
end

close all;

figure();
subplot(2,1,1)
plot(energy_kmeans_array)
title('energy')
subplot(2,1,2)
plot(condition_number_kmeans_array)
title('condition number')


figure();
for i = 1 : 3
    subplot(3,2,2*i-1)
    hist(masses_kmeans(:,1,i), 20);
    xmin = min(masses_kmeans(:,1,i));
    xmax = max(masses_kmeans(:,1,i));
    xlim([xmin, xmax]);
    title('iteration 1');
    subplot(3,2,2*i);
    hist(masses_kmeans(:,end,i), 20);
    title(['iteration ', num2str(n_iterations)]);
    xlim([xmin, xmax]);
end
