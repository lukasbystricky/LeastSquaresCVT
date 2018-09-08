function in = total_degree_set(n)

in_tmp = zeros((n+1)^2,2);
degrees = 0:n;

for i = 1:n+1
    in_tmp((i-1)*(n+1) + 1: i*(n+1),1) = degrees(i);
    for j = 1:n+1
        in_tmp((i-1)*(n+1) + j ,2) = degrees(j);
    end
end

in = [];
for i = 1:(n+1)^2
   if  sum(in_tmp(i,:)) <= n
       in(end+1,:) = in_tmp(i,:);
   end
end

