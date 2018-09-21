function in = total_degree_set(n, d)

in = zeros(1, d);

for i = 1 : n
    in(end + 1, :) = [zeros(1, d - 1), i];
    
    while norm(in(end, :) - [i, zeros(1, d - 1)]) > 0
        in(end + 1, :) = mono_total_next_grlex(d, i, in(end,:));
    end
end

