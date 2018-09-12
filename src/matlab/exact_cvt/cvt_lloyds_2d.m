function [generators, masses, energy] = cvt_lloyds_2d(generators0, n_iter, f)

n = size(generators0, 1);
generators = zeros(n, 2, n_iter);
generators(:,:,1) = generators0;
masses = zeros(n,n_iter);
energy = zeros(1, n_iter);


for i = 1:n_iter
    
    %% Update generators
    [generators(:,:,i + 1), masses(:,i), energy(i)] = ...
                    cvt_step_lloyds_2d(generators(:,:,i), f);
                
    mu = mean(masses(:, i));
    sigma = std(masses(:, i));
    
    %% Plot Voronoi diagram and mass histogram on previous diagram
    figure(1);
    clf
    
    voronoi(generators(:,1,i),generators(:,2,i));
    title(['iteration ', num2str(i)]);
    
    hold on
    for g = 1:n
        r = 0.025;
        rg = r*(masses(g, i)/mu);
        xs = generators(g, 1, i) - r;
        ys = generators(g, 2, i) - r;
        
        xsg = generators(g, 1, i) - rg;
        ysg = generators(g, 2, i) - rg;
        
        if masses(g, i) > mu
            color = 'b';
        else
            color = 'g';
        end
        
        rectangle('Position', [xsg, ysg, 2*rg, 2*rg], 'Curvature', [1,1], 'FaceColor', color);
        rectangle('Position', [xs, ys, 2*r, 2*r], 'Curvature', [1,1], 'EdgeColor', 'r');
    end
    
    axis equal
    xlim([-1,1]);
    ylim([-1,1]);

    
    figure(2);
    
    hist(masses(:,i), 20);
    title(['iteration ', num2str(i), ', mu = ', num2str(mu), ', sigma = ', num2str(sigma)]);    
    drawnow
   
end
