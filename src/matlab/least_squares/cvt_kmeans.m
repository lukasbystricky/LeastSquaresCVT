function [gens] = cvt_kmeans(gens, Phi, num_sample)

[n,~] = size(gens);
[m,d] = size(Phi.index_set);

Pweight = @(x) sum(Phi.value(x).^2, 2) / m;

scale = (d+2.0)/(d*1.0);

sample = sequential_sampling_uniform(num_sample,Phi);
weights = zeros(num_sample, 1);

for i=1:num_sample
    weights(i) = ( Pweight(sample(i,:))/(2^d)).^(scale - 1);
end
    
iterations = 100;

for iter=1:iterations
    energy = 0; 
    bins = zeros(n, d);
    bin_count = zeros(n, 1);
    
    for i=1:num_sample
       [~, ind] = min( sum((gens - sample(i,:)).^2, 2) );
       bins(ind,:) = bins(ind,:) + sample(i,:)*weights(i);
       bin_count(ind) = bin_count(ind) + weights(i);
    end
    
    for j=1:n
       if bin_count(j) ~= 0
           gens(j,:) = bins(j,:) ./ bin_count(j);
       end
    end
end

end