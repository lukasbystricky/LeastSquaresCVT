function generators0 = initial_generators(n, M, f)

generators0 = zeros(n, 2);

for i = 1 : n
    
    x = 2 * rand(1,1) - 1;
    y = 2 * rand(1,1) - 1;
    
    mu = M * rand(1,1);
    
    while f(x,y) < mu
        x = 2 * rand(1,1) - 1;
        y = 2 * rand(1,1) - 1;
        
        mu = M * rand(1,1);
    end
    
    generators0(i,:) = [x, y];
end