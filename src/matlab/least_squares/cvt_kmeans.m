function [generators,energy, masses] = cvt_kmeans(generators0, Phi, ...
                    num_sample, sample_type, iterations, adaptive, save_name)
%% cvt_kmeans constructs an approximate CVT of the Christoffel function 
%
% Creates an approximation of to a CVT of the Christoffel funtion on the
% hypercube [-1,1]^d. Instead of constructing an actual CVT, instead a
% number of possibly weighted points are placed randomly inside the domain 
% and k-means is used to find the optimal generators. 
%
% PARAMETERS:
% INPUTS:
% generators0, REAL (n,d), the initial generators
% Phi, struct containing the following information:
%   - index set, INTEGER (m,d), the index set needed to construct the basis
%   - basis_card, INTEGER, the cardinality of the index set (i.e. m)
%   - value, FUNCTION HANDLE, a function handle that specifies how to
%   compute the basis, for example @(X) eval_phi_fast(X, Phi.index_set)
% num_sample, INTEGER, the number of sample points placed in the domain
% sample_type, STRING, specifies how to draw the samples. Must be one of:
%   - 'sequential', draw the sample from the Christoffel function using the
%   sequential sampling technique described in "Optimal weighted 
%   least-squares methods" by Cohen and Migliorati (2016)
%  - 'cvt', draw the samples from the corrected CVT density, i.e.
%  \rho^(d+2/d), where \rho is the Christoffel function. This is done using
%  rejection sampling and is quite slow in higher dimensions
%  - 'uniform', draw the points from the uniform distribution on [-1,1]^d
% iterations, INTEGER, number of k-means iterations
% adaptive, BOOLEAN, if TRUE compute the change in energy each iteration,
%   if it is below a certain tolerance, double the number of sample points
%   for the next iteration
%
% OUTPUT:
% generators, REAL (n,d,iterations), the k-means generators at each
%   iteration
% energy, REAL (iterations), the k-means energy at each iteration
% masses, REAL (n, iterations), the mass of each k-means cluster at each
%  iteration

[n_generators, ~] = size(generators0);
[m, d] = size(Phi.index_set);

generators = zeros(n_generators, d, iterations);
generators(:, :, 1) = generators0;

cvt_scale = (d + 2.0)/(d * 1.0);
Pweight = @(x) sum(Phi.value(x).^2, 2) / m;

weights = zeros(num_sample, 1);

switch sample_type
    case 'sequential'
        
        % TO FIX, THIS FUNCTION IS FAR TOO SLOW
        sample = sequential_sampling_uniform(num_sample,Phi);
        
        for i = 1 : num_sample
            weights(i) = ( Pweight(sample(i,:))/(2^d)).^(cvt_scale - 1);
        end
        
    case 'cvt'
        
        sample = rejection_sampling_uniform(num_sample,Phi, true);        
        weights = ones(num_sample,1);
        
    case 'uniform'
        
        sample = zeros(num_sample, d);
        for k = 1 : d
            sample(:, k) = 2 * rand(num_sample, 1) - 1;
        end
        
        for i = 1 : num_sample
           weights(i) = Pweight(sample(i,:)).^cvt_scale;
        end
end

energy = zeros(iterations, 1);
masses = zeros(n_generators, iterations);

for iter = 1 : iterations
    
        
    if d == 2
        [~, ~, energy(iter)]  =  vornoi_compute_mass(generators(:,:,iter), ...
            @(x) (sum(Phi.value(x).^2, 2) / m / (2^d)));
    else
        energy(iter) = 0;
    end

   [k, dist] = dsearchn(generators(:,:,iter), sample);
    bins = zeros(n_generators, d);
    bin_count = zeros(n_generators, 1);
    
    for i = 1 : n_generators
        bins(i,:) = sum(bsxfun(@times, sample(k == i,:), weights(k == i)),1);
        bin_count(i) = sum(weights(k == i));
        
        masses(i, iter) = sum(Pweight(sample(k == i,:)))/ num_sample;
        
        if d ~= 2
            energy(iter) = energy(iter) + sum(weights(k == i) .* ...
                dist(k == i).^2) / num_sample;
        end
    end
  
%     for i = 1 : num_sample
%         [~, ind] = min( sum((generators(:, :, iter)...
%             - repmat(sample(i, :), n_generators, 1)).^2, 2) );
%         bins(ind,:) = bins(ind,:) + sample(i, :) * weights(i);
%         bin_count(ind) = bin_count(ind) + weights(i);
%         
%         masses(ind, iter) = masses(ind, iter) + Pweight(sample(i, :)) ...
%             / num_sample;
%         
%         if d ~= 2
%             energy(iter) = energy(iter) + weights(i) ...
%                 * norm(generators(ind, :, iter) - sample(i, :))^2 ...
%                 / num_sample;
%         end
%     end
    
    disp( ['iteration ',num2str(iter), ' energy: ', num2str(energy(iter))] )
    
    for j = 1 : n_generators
       if bin_count(j) ~= 0
           generators(j,:,iter + 1) = bins(j,:) ./ bin_count(j);
       end
    end
    
    if ~isempty(save_name)
       save(save_name, 'generators', 'energy'); 
    end
    
    % if energy change is small, double number of sample points
    if adaptive && iter > 5 ...
            && abs((energy(iter) - energy(iter - 1))/energy(iter)) < 1e-4 && num_sample < 5e6
        
        disp(['Doubling sample points to ', num2str(2* num_sample)]);
        
        sample = [sample; zeros(num_sample, d)];
        weights = [weights; zeros(num_sample, 1)];
        
        switch sample_type
            case 'sequential'
                
                % TO FIX, THIS FUNCTION IS FAR TOO SLOW
                sample(num_sample + 1 : end) ...
                            = sequential_sampling_uniform(num_sample,Phi);
                
                for i = num_sample + 1 : 2 * num_sample
                    weights(i) = ( Pweight(sample(i,:)) ...
                                / (2^d)).^(cvt_scale - 1);
                end
                
            case 'cvt'
                
                % TO IMPLEMENT
                
            case 'uniform'
                
                for k = 1 : d
                    sample(num_sample + 1 : end, k) = 2 * ...
                                rand(num_sample, 1) - 1;
                end
                
                for i = num_sample + 1 : 2 * num_sample
                    weights(i) = Pweight(sample(i,:)).^cvt_scale;
                end
        end
        
        num_sample = 2* num_sample;
        
        
    end
end