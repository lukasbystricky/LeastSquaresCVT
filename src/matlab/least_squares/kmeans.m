function [gens, prevWeights] = kmeans(gens, Phi, num_sample, prevWeights)

prevGens = gens;
prevEnergy = Inf;

[n,~] = size(gens);
[m,d] = size(Phi.index_set);

Pweight = @(x) sum(Phi.value(x).^2, 2) / m;

scale = (d+2.0)/(d*1.0);

sample = rejection_sampling_uniform(num_sample,Phi);
weights = zeros(num_sample, 1);

for i=1:num_sample
    weights(i) = ( Pweight(sample(i,:))/(2^d)).^(scale - 1);
end
    
for k=1:25
    energy = 0; 
    bins = zeros(n, d);
    bin_count = zeros(n, 1);
    
    for i=1:num_sample
       [~, ind] = min( sum((gens - sample(i,:)).^2, 2) );
       bins(ind,:) = bins(ind,:) + sample(i,:)*weights(i);
       bin_count(ind) = bin_count(ind) + weights(i);
       energy = energy + weights(i)*norm( sample(i,:) - gens(ind) );
    end
    
    if abs(prevEnergy - energy) < 1e-8
        break
    end
    
    prevEnergy = energy;
    
    for j=1:n
       if bin_count(j) ~= 0
           gens(j,:) = ( bins(j,:) + prevGens(j,:)*prevWeights(j) ) ...
                        ./ ( bin_count(j) + prevWeights(j));
       end
    end
    
    prevWeights = prevWeights + bin_count;
    prevGens = gens;
end

end