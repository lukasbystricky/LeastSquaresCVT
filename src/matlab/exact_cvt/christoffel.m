function W = christoffel(total_degree, X)

x = X(:,1);
y = X(:,2);

index_set = total_degree_set(total_degree);

m = size(index_set,1);
drho = 0.25;

W = 0;
for j = 1:m
   W_tmp = sqrt(1 + 2*index_set(j,1))* sqrt(1 + 2*index_set(j,2)) * ...
        polyval(LegendrePoly(index_set(j,1)),x) .* polyval(LegendrePoly(index_set(j,2)),y) ;
   W = W + W_tmp.^2 * drho;
end

W = W/m;