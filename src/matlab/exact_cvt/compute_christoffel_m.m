function M = compute_christoffel_m(td)

index_set = total_degree_set(td);

m = size(index_set, 1);

M = 0;

for i = 1 : m
    M = M + (1 + 2*index_set(i,1))*(1 + 2*index_set(i,2));
end
