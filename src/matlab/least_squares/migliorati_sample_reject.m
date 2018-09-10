function omega = migliorati_sample_reject(Omega, weight, Napprox)

% generate N samples from density rho*weight,
% where rho is the uniform density

xl = Omega.box(1,:);
xr = Omega.box(2,:);

% the weight is maximal in a corner of the cube
% this is just a conjecture
w_max = weight(xl);

Ntrial = ceil(Napprox * w_max);

Xtrial = bsxfun(@plus,xl, (bsxfun(@times, (xr-xl), rand(Ntrial, Omega.d))));
wXtrial = weight(Xtrial);

assert(all(wXtrial <= w_max))

accept = w_max*rand(Ntrial,1) < wXtrial;

omega = struct();
omega.x = Xtrial(accept,:);
omega.lambda = 1 ./ (wXtrial(accept) * length(find(accept)));

end

