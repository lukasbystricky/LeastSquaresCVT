function [sample_points] = rejection_sampling_uniform(n,Phi,scaled)
%% rejection_sampling_uniform creates sample points with rejection sampling
%
% Creates sample points according to the Christoffel function, possibly
% raised to the power of (d+2)/d, the appropriate CVT scaling factor. Does
% so with rejection sampling, in which samples from max{d\rho} are sampled,
% then rejected if the value is greater than d\rho.
%
% PARAMETERS:
% INPUTS:
% n, INTEGER, the number of sample poitns to produce
% Phi, struct containing the following information:
%   - index set, INTEGER (m,d), the index set needed to construct the basis
%   - basis_card, INTEGER, the cardinality of the index set (i.e. m)
%   - value, FUNCTION HANDLE, a function handle that specifies how to
%   compute the basis, for example @(X) eval_phi_fast(X, Phi.index_set)
% scaled, BOOLEAN, if true, evaluate density raised to (d+2)/d
% OUTPUT:
% sample_points, real(n,d), the resulting sample points

[m,d] = size(Phi.index_set);

Pweight = @(x) sum(Phi.value(x).^2, 2);
max_p = Pweight(ones(1,d))/m/2^d;

failed = ones(15);
counter = 1;
sample_points = zeros(n, d);

while counter ~= n + 1
    ys = rand(1000,d).*2 - 1;
    tests = rand(1000,1).*max_p;

    for i=1:1000
       if ( Pweight(ys(i,:))/m/2^d ) >= tests(i) && counter ~= n + 1
           sample_points(counter,:) = ys(i,:);
           counter = counter + 1;
       else
           failed = failed + 1;
       end
    end
end

end