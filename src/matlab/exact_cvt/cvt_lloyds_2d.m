function [generators, masses] = cvt_lloyds_2d(generators0, n_iter, f)

n = size(generators0, 1);
generators = zeros(n, 2, n_iter);
generators(:,:,1) = generators0;
masses = zeros(n,n_iter);

for i = 1:n_iter
    [generators(:,:,i + 1), masses(:,i + 1)] = cvt_step_lloyds_2d(generators(:,:,i), f);
    mu = mean(masses(:,i + 1));
    sigma = std(masses(:,i + 1));
    
    
    figure(1);
    clf
    
    voronoi(generators(:,1,i+1),generators(:,2,i+1));
    title(['iteration ', num2str(i)]);
    
    hold on
    for g = 1:n
        r = 0.025;
        rg = r*(masses(g, i + 1)/mu);
        xs = generators(g,1,i+1) - r;
        ys = generators(g,2,i+1) - r;
        
        xsg = generators(g,1,i+1) - rg;
        ysg = generators(g,2,i+1) - rg;
        
        if masses(g, i + 1) > mu
            color = 'b';
        else
            color = 'g';
        end
        
        rectangle('Position', [xsg, ysg, 2*rg, 2*rg], 'Curvature', [1,1], 'FaceColor', color);
        rectangle('Position', [xs, ys, 2*r, 2*r], 'Curvature', [1,1], 'EdgeColor', 'r');
    end
    
    xlim([-1,1]);
    axis equal

    
    figure(2);
    
    hist(masses(:,i + 1), 20);
    title(['iteration ', num2str(i), ', mu = ', num2str(mu), ', sigma = ', num2str(sigma)]);
   % xlim([0,1]);
    
    drawnow
end
