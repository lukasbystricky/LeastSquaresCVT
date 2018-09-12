function cvt_animation(generators, masses, fname)

n_iter = size(generators,3);
n_generators = size(generators,1);

h = figure();

for i = 1:n_iter
    clf
    voronoi(generators(:,1,i), generators(:,2,i));
    hold on
    
    mu = mean(masses(:,i));
    
    for g = 1:n_generators
        r = 0.025;
        rg = r*(masses(g, i)/mu);
        xs = generators(g,1,i) - r;
        ys = generators(g,2,i) - r;
        
        xsg = generators(g,1,i) - rg;
        ysg = generators(g,2,i) - rg;
        
        if masses(g, i) > mu
            color = 'b';
        else
            color = 'g';
        end
        
        rectangle('Position', [xsg, ysg, 2*rg, 2*rg], 'Curvature', [1,1],...
                            'FaceColor', color);
        rectangle('Position', [xs, ys, 2*r, 2*r], 'Curvature', [1,1],...
                            'EdgeColor', 'r');
    end
    
    title(['iteration ', num2str(i)]);
    
    axis equal
    xlim([-1,1]);
    ylim([-1,1]);
    
    drawnow
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 2
        imwrite(imind,cm, [fname, '.gif'] ,'gif', 'Loopcount',inf, ...
                        'DelayTime', 0);
    else
        imwrite(imind,cm, [fname, '.gif'],'gif', 'WriteMode','append', ...
                        'DelayTime', 0);
    end
end
